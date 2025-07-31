import matlab.engine
from flask_cors import CORS
import os
import json
import threading
import time
from flask import Flask, send_from_directory, jsonify, request
from flask_socketio import SocketIO, emit 
import io
import base64
from PIL import Image
import numpy as np


# --- Flask and SocketIO Setup ---
app = Flask(__name__, static_folder=os.path.join(os.path.dirname(__file__), '../frontend/build'))
CORS(app)
app.config['SECRET_KEY'] = 'your_super_secret_key_here' # IMPORTANT: Change this to a strong, unique key
socketio = SocketIO(app, cors_allowed_origins="*")


# --- Global MATLAB Engine and Lock ---
eng = None # Renamed from 'eng' to 'matlab_eng' for clarity in new functions
matlab_ready_event = threading.Event()
matlab_lock = threading.Lock() # Define the lock globally

# Directory where MATLAB will save temporary batch files
MATLAB_TEMP_DIR = os.path.expanduser('~/Desktop/matlab/temp_results') # Must match MATLAB's temp_dir
MATLAB_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tempfiles')
MATLAB_OUTPUT_FILENAME = 'radar_results.json' # Updated filename
MATLAB_OUTPUT_PATH = os.path.join(MATLAB_OUTPUT_DIR, MATLAB_OUTPUT_FILENAME)
os.makedirs(MATLAB_OUTPUT_DIR, exist_ok=True)
# --- MATLAB Engine Initialization ---
def start_matlab_engine():
    global eng # Use 'eng' to match your existing code's usage
    print("Starting MATLAB engine...")
    try:
        eng = matlab.engine.start_matlab()
        matlab_files_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../matlabfiles'))
        eng.addpath(matlab_files_path, nargout=0)
        print(f"MATLAB engine started. Added path: {matlab_files_path}")
        matlab_ready_event.set()
    except Exception as e:
        print(f"Failed to start MATLAB engine: {e}")

# Start MATLAB engine in a separate thread to not block Flask startup
threading.Thread(target=start_matlab_engine, daemon=True).start()

# --- Flask Routes ---
@app.route('/')
def index():
    return send_from_directory(os.path.join(os.path.dirname(__file__), '../frontend/build'), 'index.html')

@app.route('/api/get_latest_radar_data', methods=['GET'])
def get_latest_radar_data():
    try:
        if not os.path.exists(MATLAB_OUTPUT_PATH):
            print(f"Error: JSON file not found at {MATLAB_OUTPUT_PATH}") # Add this log
            return jsonify({"error": "No radar analysis data file found. Please run an analysis first."}), 404

        with open(MATLAB_OUTPUT_PATH, 'r') as f:
            data = json.load(f)
        return jsonify(data)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {MATLAB_OUTPUT_PATH}") # Add this log
        return jsonify({"error": "Invalid JSON format in radar results file."}), 500
    except Exception as e:
        print(f"Unexpected error reading file {MATLAB_OUTPUT_PATH}: {str(e)}") # Add this log
        return jsonify({"error": f"Failed to read radar data: {str(e)}"}), 500


@app.route('/list_matlab_files')
def list_matlab_files():
    matlab_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../matlabfiles'))
    try:
        files = [f for f in os.listdir(matlab_dir) if f.endswith('.m')]
        return jsonify(files)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# --- Helper Function for MATLAB File Info ---
def get_function_info(func_name):
    """
    Determines if a MATLAB .m file is a function, a script, or neither.
    Returns a tuple: (nargout_count_if_function, is_script_flag)
    - nargout_count_if_function: The number of outputs if it's a true function, else 0.
    - is_script_flag: True if it's a script (identified as 3 by exist, or as a script by nargout error), else False.
    """
    try:
        if eng is None or not matlab_ready_event.is_set():
            print("WARNING: MATLAB engine not ready in get_function_info.")
            return (0, False)

        matlab_file_type = eng.eval(f"exist('{func_name}', 'file')", nargout=1)
        print(f"DEBUG (get_function_info): exist('{func_name}', 'file') returned {matlab_file_type}")

        if matlab_file_type == 3: # It's definitely a script (.m file)
            print(f"INFO: '{func_name}' is identified as a SCRIPT by 'exist'.")
            return (0, True) # 0 outputs, it's a script
        elif matlab_file_type == 2: # It's identified as a function (could still be a script that looks like a function)
            try:
                n_out = eng.eval(f"nargout('{func_name}')", nargout=1)
                print(f"INFO: '{func_name}' is confirmed as a FUNCTION (nargout: {n_out}).")
                return (int(max(0, n_out)), False) # True function
            except matlab.engine.MatlabExecutionError as e:
                if "is a script" in str(e):
                    print(f"INFO: nargout for '{func_name}' failed, indicating it's a SCRIPT despite 'exist' returning 2.0.")
                    return (0, True) # It's actually a script
                elif "Illegal use of reserved keyword" in str(e):
                    print(f"ERROR: '{func_name}' is malformed (Illegal use of reserved keyword): {str(e)}")
                    return (0, True) # Treat as script (best effort, expects file fixed)
                else:
                    print(f"WARNING: Unknown nargout error for '{func_name}': {str(e)}")
                    return (0, False) # Assume 0-output function to be safe
        else: # Not found or other type
            print(f"WARNING: '{func_name}' not identified as function/script (type: {matlab_file_type}). Assuming script for eval.")
            return (0, True) # Default to script to try running it with 'run()'

    except Exception as e:
        print(f"Error determining file type for {func_name}: {e}")
        return (0, False) # Default to 0 outputs, not a special script case

# --- NEW SocketIO Event Handler for running user-provided MATLAB code with variables ---
@socketio.on('run_editor_matlab_code')
def handle_run_editor_matlab_code(data):
    original_code = data.get('code', '') 
    variables = data.get('variables', {}) 
                                        
    print(f"\n>>> SOCKET EVENT: 'run_editor_matlab_code' received.")
    print(f"Original Code Length: {len(original_code)} characters")

    if not matlab_ready_event.wait(timeout=30):
        error_msg = 'MATLAB engine not ready after timeout. Please restart server.'
        print(error_msg)
        socketio.emit('matlab_code_output', {'error': error_msg})
        return
    if eng is None:
        error_msg = 'MATLAB engine is not initialized.'
        print(error_msg)
        socketio.emit('matlab_code_output', {'error': error_msg})
        return

    with matlab_lock:
        stdout_stream = io.StringIO()
        stderr_stream = io.StringIO()
        try:
        
            print(f"Executing MATLAB code from editor (variables are embedded in the code string):\n--- START EDITOR CODE ---\n{original_code}\n--- END EDITOR CODE ---")
            eng.eval(original_code, nargout=0, stdout=stdout_stream, stderr=stderr_stream)

            output_lines = stdout_stream.getvalue().strip()
            error_lines = stderr_stream.getvalue().strip()

            response_data = {}
            if output_lines:
                response_data['output'] = output_lines.split('\n')
                print(f"MATLAB Standard Output:\n{output_lines}")
            if error_lines:
                response_data['error'] = error_lines
                print(f"MATLAB Standard Error:\n{error_lines}")

            if not response_data:
                response_data['output'] = ["MATLAB script executed with no output."]

            socketio.emit('matlab_code_output', response_data)

        except matlab.engine.MatlabExecutionError as me:
            print(f"MATLAB execution error: {me.args[0]}")
            socketio.emit('matlab_code_output', {'error': f"MATLAB Execution Error: {me.args[0]}"})
        except Exception as e:
            print(f"An unexpected error occurred during MATLAB execution: {e}")
            socketio.emit('matlab_code_output', {'error': f"Server Error: {str(e)}"})

# --- Existing SocketIO Event Handler for running existing .m files (e.g., from Dashboard) ---
@socketio.on('run_matlab_code')
def run_matlab_code(data):
    print("\n>>> SOCKET EVENT: 'run_matlab_code' received (for existing .m files).")
    global eng
    if not matlab_ready_event.wait(timeout=30):
        socketio.emit('matlab_result', {'error': 'MATLAB engine not ready after timeout. Please restart server.'})
        return
    if eng is None:
        socketio.emit('matlab_result', {'error': 'MATLAB engine is not initialized.'})
        return

    filename = data.get('filename')
    args = data.get('args', ()) # This 'args' contains the dictionary of parameters

    if not filename:
        socketio.emit('matlab_result', {'error': 'No MATLAB file specified.'})
        return

    func_name = os.path.splitext(filename)[0]
    print(f"Attempting to run existing MATLAB file: {filename}")
    print(f"Extracted func_name: '{func_name}'")
    print(f"Arguments for MATLAB function: {args}")

    try:
        results_list = []
        nargout_count, is_really_script = get_function_info(func_name)

        with matlab_lock:
            stdout_stream = io.StringIO()
            stderr_stream = io.StringIO()

            if is_really_script:
                print(f"'{func_name}' is being handled as a SCRIPT. Running it using eng.eval()...")
                # MODIFIED START: Correctly access the parameters from the 'args' list
                if args and isinstance(args, (list, tuple)) and len(args) > 0 and isinstance(args[0], dict):
                    config_data = args[0] # Get the dictionary of parameters
                    print(f"Injecting configuration variables into MATLAB workspace: {config_data}")
                    for var_name, var_value in config_data.items():
                        try:
                            # Attempt to convert string values to appropriate Python types
                            if isinstance(var_value, str):
                                if var_value.lower() == 'true':
                                    python_value = True
                                elif var_value.lower() == 'false':
                                    python_value = False
                                else:
                                    try:
                                        # Try converting to float, then int if it's a whole number
                                        float_value = float(var_value)
                                        if float_value.is_integer():
                                            python_value = int(float_value)
                                        else:
                                            python_value = float_value
                                    except ValueError:
                                        # If not a boolean or number, keep as string
                                        python_value = var_value
                            else:
                                python_value = var_value # Use value as is if not a string

                            eng.workspace[var_name] = python_value
                            print(f"    - Set '{var_name}' = '{python_value}' (Python type: {type(python_value)}) in MATLAB workspace.")
                        except Exception as e:
                            print(f"    - ERROR injecting variable '{var_name}' with value '{var_value}': {e}")
                else:
                    print("No configuration data (parameters) found in 'args' to inject for the script.")
 

                eng.eval(f"run('{filename}')", nargout=0, stdout=stdout_stream, stderr=stderr_stream);
 
        output_lines_from_file = stdout_stream.getvalue().strip()
        error_lines_from_file = stderr_stream.getvalue().strip()

        results_converted = [convert_output(r) for r in results_list]

        socketio.emit('matlab_result', {
            'results': results_converted,
            'message': f'Ran {func_name} successfully with {len(results_converted)} outputs.',
            'stdout': output_lines_from_file.split('\n') if output_lines_from_file else [],
            'stderr': error_lines_from_file
        })

    except matlab.engine.MatlabExecutionError as e:
        print(f"MATLAB Execution Error in run_matlab_code for {filename}: {str(e)}")
        socketio.emit('matlab_result', {'error': str(e)})
    except Exception as e:
        print(f"Server error in run_matlab_code for {filename}: {type(e).__name__}: {str(e)}")
        socketio.emit('matlab_result', {'error': str(e)})

# --- MODIFIED: This section will now handle the full simulation including batching ---
@socketio.on('run_matlab_simulation') # Frontend will trigger this event
def handle_run_matlab_simulation(data):
    global eng # Use 'eng' to match your existing code's usage
    print("Received request to run MATLAB simulation.")

    # Optional: Clear temp results directory before starting a new run
    if os.path.exists(MATLAB_TEMP_DIR):
        for f in os.listdir(MATLAB_TEMP_DIR):
            os.remove(os.path.join(MATLAB_TEMP_DIR, f))
        print(f"Cleaned up {MATLAB_TEMP_DIR}")
    else:
        os.makedirs(MATLAB_TEMP_DIR)
        print(f"Created {MATLAB_TEMP_DIR}")

    if eng is None: # Check if engine is already running from startup thread
        try:
            # This block should ideally not be hit if start_matlab_engine() runs on startup
            # But it's a fallback if the startup thread fails or engine quits.
            eng = matlab.engine.start_matlab()
            print("MATLAB engine started (fallback).")
            matlab_files_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../matlabfiles'))
            eng.addpath(matlab_files_path, nargout=0)
            print(f"MATLAB engine started. Added path: {matlab_files_path}")
            matlab_ready_event.set()
        except Exception as e:
            print(f"Failed to start MATLAB engine: {e}")
            emit('simulation_error', {'message': f'Failed to start MATLAB engine: {e}'})
            return
    
    # Ensure engine is ready before proceeding
    if not matlab_ready_event.wait(timeout=30):
        emit('simulation_error', {'message': 'MATLAB engine not ready after timeout. Please restart server.'})
        return


    with matlab_lock:
        stdout_stream = io.StringIO()
        stderr_stream = io.StringIO()

        # Start a thread to read MATLAB output as it comes
        # This thread will also trigger process_matlab_batch_file
        output_reader_thread = threading.Thread(target=read_matlab_output, args=(stdout_stream, socketio))
        output_reader_thread.daemon = True # Thread will exit when main program exits
        output_reader_thread.start()

        try:
            matlab_script_name = 'unified_usrp_test_30k.m'
            matlab_script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../matlabfiles', matlab_script_name))

            # Set MATLAB's current directory to where the script is located
            eng.cd(os.path.dirname(matlab_script_path))
            print(f"MATLAB current directory set to: {os.path.dirname(matlab_script_path)}")

            # Inject configuration variables from the frontend (if any)
            # The 'data' payload for this event is the config, e.g., {'runTests': True, 'sensor_test': False, ...}
            print(f"Injecting configuration variables from frontend into MATLAB workspace: {data}")
            for var_name, var_value in data.items():
                try:
                    if isinstance(var_value, str) and var_value.lower() in ('true', 'false'):
                        python_value = var_value.lower() == 'true'
                    elif isinstance(var_value, str):
                        try:
                            python_value = float(var_value)
                            if python_value.is_integer():
                                python_value = int(python_value)
                        except ValueError:
                            python_value = var_value
                    else:
                        python_value = var_value

                    eng.workspace[var_name] = python_value
                    print(f"    - Set '{var_name}' = '{python_value}' (Python type: {type(python_value)}) in MATLAB workspace.")
                except Exception as e:
                    print(f"    - ERROR injecting variable '{var_name}' with value '{var_value}': {e}")

            print(f"Calling MATLAB script: {matlab_script_name}")
            eng.eval(f"run('{matlab_script_name}')", nargout=0, stdout=stdout_stream, stderr=stderr_stream)
            print("MATLAB script finished execution.")

            # After the MATLAB script completes, you can optionally send a 'finished' message
            # The read_matlab_output thread will eventually finish as well when stdout closes
            socketio.emit('simulation_complete', {'message': 'MATLAB simulation finished. All batches processed.'})

        except Exception as e:
            print(f"Error running MATLAB script: {e}")
            emit('simulation_error', {'message': f'Error during MATLAB simulation: {e}'})
        finally:
            # Ensure streams are closed so the reading thread can terminate
            stdout_stream.close()
            stderr_stream.close()
            # The output_reader_thread is a daemon thread, so it will exit with the main program.
            # If it's not daemon, you'd need output_reader_thread.join() here.

# --- Original 'run_radar_analysis' is now commented out as 'run_matlab_simulation' handles the full workflow ---
# @socketio.on('run_radar_analysis')
# def handle_run_radar_analysis(data):
#     # ... (original code for run_radar_analysis) ...
#     pass # This event handler is now superseded by 'run_matlab_simulation' for unified_usrp_test_30k.m

# --- Core Helper for Converting MATLAB Types to Python Types ---
def convert_output(output):
    """
    Recursively converts MATLAB engine outputs to standard Python types.
    Handles native Python types first for robustness.
    """
    # 1. Handle native Python types (float, int, bool, str, list, tuple, dict)
    if isinstance(output, (float, int, bool, str)):
        return output
    elif isinstance(output, (list, tuple)):
        return [convert_output(item) for item in output]
    elif isinstance(output, dict):
        py_dict = {}
        for key, value in output.items():
            py_dict[key] = convert_output(value)
        return py_dict

    # 2. Then handle matlab.engine specific types
    elif isinstance(output, matlab.double):
        # print(f"DEBUG: convert_output received matlab.double object. Value: {output}, Raw type: {type(output)}")

        try:
            if hasattr(output, 'tolist') and callable(getattr(output, 'tolist')):
                python_list = output.tolist() # Get the Python list representation

                # Apply flattening logic here for common MATLAB vector formats
                if isinstance(python_list, list) and all(isinstance(sublist, list) and len(sublist) == 1 for sublist in python_list):
                    # Column vector like [[val],[val]] becomes [val,val]
                    # print("DEBUG: Flattening list of single-element lists (e.g., [[val], [val]]).")
                    return [item[0] for item in python_list]
                elif isinstance(python_list, list) and len(python_list) == 1 and isinstance(python_list[0], list):
                    # Row vector like [[val1, val2, ...]] becomes [val1, val2, ...]
                    # print("DEBUG: Flattening single-list-wrapped list (e.g., [[val1, val2]]).")
                    return python_list[0]
                else:
                    # Multi-dimensional arrays or already flat lists
                    # print("DEBUG: matlab.double converted to Python list (no further flattening needed).")
                    return python_list
            else: # tolist is missing or not callable, proceed to manual iteration fallback
                print(f"CRITICAL FALLBACK: matlab.double (type: {type(output)}) has NO .tolist() or .size. Attempting manual iteration and float conversion of elements.")
                converted_data = []
                try:
                    temp_elements = []
                    try:
                        for item in output:
                            temp_elements.append(item)
                    except TypeError:
                        try:
                            print("DEBUG: matlab.double not iterable, trying direct float conversion as scalar.")
                            return float(output)
                        except (TypeError, ValueError):
                            print(f"CRITICAL ERROR: matlab.double object is neither iterable nor directly float-convertible, and missing size/tolist. Cannot convert. Object: {output}")
                            return output

                    if not temp_elements:
                        print("DEBUG: matlab.double iterable but empty.")
                        return []

                    if isinstance(temp_elements[0], matlab.double):
                        print("DEBUG: matlab.double is a nested array (matrix). Processing rows.")
                        for row_obj in temp_elements:
                            sub_list = []
                            for val in row_obj:
                                try:
                                    sub_list.append(float(val))
                                except (TypeError, ValueError):
                                    print(f"ERROR: Cannot convert sub-element '{val}' (type: {type(val)}) to float during nested iteration. Returning raw.")
                                    sub_list.append(val)
                            converted_data.append(sub_list)
                    else:
                        print("DEBUG: matlab.double is a flat array (vector or scalar). Processing elements.")
                        for val in temp_elements:
                            try:
                                converted_data.append(float(val))
                            except (TypeError, ValueError):
                                print(f"ERROR: Cannot convert element '{val}' (type: {type(val)}) to float during flat iteration. Returning raw.")
                                converted_data.append(val)
                    return converted_data

                except Exception as e:
                    print(f"CRITICAL ERROR during manual matlab.double iteration/conversion (inner catch): {e}. Object: {output}. Returning raw.")
                    return output

        except Exception as e:
            print(f"CRITICAL ERROR in matlab.double conversion (outer catch for standard attributes): {e}. Object type: {type(output)}. Returning as-is.")
            return output

    elif isinstance(output, matlab.logical):
        return bool(output)

    # 3. Handle Structs and CellArrays with a protective try-except
    try:
        if isinstance(output, eng.StructArray):
            py_list = []
            for s in output.tolist():
                py_dict = {}
                for field_name in s._fieldnames:
                    py_dict[field_name] = convert_output(getattr(s, field_name))
                py_list.append(py_dict)
            return py_list
        elif isinstance(output, eng.Struct):
            py_dict = {}
            for field_name in output._fieldnames:
                py_dict[field_name] = convert_output(getattr(output, field_name))
            return py_dict
        elif isinstance(output, matlab.cell.CellArray):
            return [convert_output(c) for c in output.tolist()]
    except AttributeError as e:
        print(f"WARNING: Could not check MATLAB specific type due to AttributeError (possible matlab.engine bug): {e}")
        print(f"DEBUG: Output type causing this specific AttributeError: {type(output)}")
        return output
    except Exception as e:
        print(f"WARNING: Unexpected error during MATLAB type check: {type(e).__name__}: {e}")
        print(f"DEBUG: Output type causing this unexpected error: {type(output)}")
        return output

    # 4. Fallback for truly unhandled types
    print(f"DEBUG: convert_output received unhandled type (returning as-is): {type(output)}")
    return output

# --- Main Application Run Block ---
if __name__ == '__main__':
    print("Starting Flask-SocketIO server...")
    socketio.run(app, debug=False, port=5000, allow_unsafe_werkzeug=True)
