import React, { useState, useEffect, useRef, useCallback } from "react";
import { io } from 'socket.io-client'; // Direct import of io
import "../App.css"; // Your existing global CSS

//For popup after matlab code is done running
import { ToastContainer, toast } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';

// CodeEditor component, now directly responsible for its MATLAB output display
// It receives setMatlabResults as a prop to update parent state if needed.
function CodeEditor({ setMatlabResults }) { // Removed setMessage prop as it's not used internally for direct output
    const [lightMode, setLightMode] = useState(true);
    const [originalCode, setOriginalCode] = useState(""); // Stores the content of the uploaded file
    const [fileName, setFileName] = useState("No file uploaded");
    const [fileContent, setFileContent] = useState(""); // The code displayed in the <pre> tag, updated with variables
    const [variables, setVariables] = useState({
        // IMPORTANT: These should ideally reflect the data types and values MATLAB expects.
        // For strings, if MATLAB expects them quoted, ensure they are.
        // For booleans, 'true'/'false' as strings often works with backend conversion.
        runTests: "true",
        sensor_test: "true",
        vital_sign_test: "true",
        speed_test: "false",
        real_time: "true",
        // Correcting the initial values based on common MATLAB usage for these names:
        ts: "0.06", // Assuming this is a numerical value (time step), should be a string here
        img_filename: "'/home/xzhang/Desktop/matlab/images/jarcat.jpg'", // Path as a string, quoted for MATLAB
        strCOM: "'/dev/ttyACM1'", // COM port as a string, quoted for MATLAB
        signal_type: "'5G'", // String value, quoted for MATLAB
        Bandwidth: "'100 MHz'", // String value, quoted for MATLAB
    });
    const [connectionStatus, setConnectionStatus] = useState("Disconnected");
    const [isLoading, setIsLoading] = useState(false); // New state for loading indicator
    const [saveStatus, setSaveStatus] = useState(""); // e.g., "Saved", "Saving..."
    const [runStatus, setRunStatus] = useState(""); // "", "success", "error"

    // Removed scanStatus, uptime, startTimeRef, pausedDurationRef from here
    // as they seemed more relevant to Dashboard and are not used in editor's core logic.
    // If you need them for general app status, they should be managed in a parent component like App.js.
    const [scanStatus, setScanStatus] = useState("Stopped"); // Keeping as per your original file
    const [uptime, setUptime] = useState("00:00:00"); // Keeping as per your original file
    const startTimeRef = useRef(null); // Keeping as per your original file
    const pausedDurationRef = useRef(0); // Keeping as per your original file


    useEffect(() => {
    const savedFileName = localStorage.getItem("uploadedFileName");
    const savedOriginalCode = localStorage.getItem("uploadedOriginalCode");
    const savedFileContent = localStorage.getItem("uploadedFileContent");

    if (savedFileName && savedOriginalCode && savedFileContent) {
        setFileName(savedFileName);
        setOriginalCode(savedOriginalCode);
        setFileContent(savedFileContent);
    }
}, []);

    // Define ref for WebSocket connection
    const socketRef = useRef(null);

    // Debounce function to avoid unnecessary regeneration
    const debounce = (func, delay) => {
        let timeout;
        return (...args) => {
            clearTimeout(timeout);
            timeout = setTimeout(() => func(...args), delay);
        };
    };

    // Debounced version of regenerateCodeWithVariables
    // This useCallback ensures the debounced function itself is stable across renders.
    const regenerateCodeWithVariablesDebounced = useCallback(
        debounce((updatedVars) => {
            regenerateCodeWithVariables(updatedVars);
        }, 500),
        [] // No dependencies, as regenerateCodeWithVariables is also memoized
    );

    // Keeping updateUptime as per your original file, though it might not be relevant here.
    const updateUptime = () => {
        if (!startTimeRef.current) return;
        let now = new Date();
        let elapsed = now - startTimeRef.current - pausedDurationRef.current;

        let hours = Math.floor(elapsed / (1000 * 60 * 60));
        let minutes = Math.floor((elapsed % (1000 * 60 * 60)) / (1000 * 60));
        let seconds = Math.floor((elapsed % (1000 * 60)) / 1000);

        setUptime(
            String(hours).padStart(2, "0") +
            ":" +
            String(minutes).padStart(2, "0") +
            ":" +
            String(seconds).padStart(2, "0")
        );
    };

    const handleFileUpload = (e) => {
        const file = e.target.files[0];
        if (file) {
            setFileName(file.name);
            const reader = new FileReader();
           reader.onload = (event) => {
    const content = event.target.result;
    setOriginalCode(content); // Store original code
    setFileContent(content); // Immediately show uploaded code (before injection)

    // Inject variables into code right after upload
    const updatedVars = { ...variables };
    let newCode = content;
    Object.entries(updatedVars).forEach(([key, val]) => {
        const regex = new RegExp(`(^|\\s)${key}\\s*=\\s*[^;\\n]*(;|$)?`, 'gm');
        let formattedVal = val;
        if (typeof val === 'string' && !val.startsWith("'") && !val.endsWith("'") && isNaN(Number(val))) {
            formattedVal = `'${val}'`;
        }
        if (newCode.match(regex)) {
            newCode = newCode.replace(regex, `$1${key} = ${formattedVal}$2`);
        }
    });
    setFileContent(newCode); // Update display

    // Save to localStorage
    localStorage.setItem("uploadedFileName", file.name);
    localStorage.setItem("uploadedOriginalCode", content);
    localStorage.setItem("uploadedFileContent", newCode);
};

            reader.readAsText(file);
        } else {
            setFileName("No file uploaded");
            setOriginalCode("");
            setFileContent("");
            // Reset variables to defaults if file is cleared
            setVariables({
                runTests: "true", sensor_test: "true", vital_sign_test: "true", speed_test: "false",
                real_time: "true", ts: "0.06", img_filename: "'/home/xzhang/Desktop/matlab/images/jarcat.jpg'",
                strCOM: "'/dev/ttyACM1'", signal_type: "'5G'", Bandwidth: "'100 MHz'",
            });
            setMatlabResults(null); // Clear MATLAB results on file clear
        }
    };

    // Memoized version of regenerateCodeWithVariables
    // This ensures the function identity is stable, which is important for useCallback dependencies.
    const regenerateCodeWithVariables = useCallback((updatedVars) => {
        let newCode = originalCode;
        if (!originalCode) { // Do not attempt to regenerate if no file is loaded
            setFileContent("");
            return;
        }

        Object.entries(updatedVars).forEach(([key, val]) => {
            // This regex attempts to find and replace existing variable assignments.
            // It looks for 'variable_name = value;' or 'variable_name = value'
            // The `g` flag ensures all occurrences are replaced.
            const regex = new RegExp(`(^|\\s)${key}\\s*=\\s*[^;\\n]*(;|$)?`, 'gm');

            let formattedVal = val;
            // Heuristic: if a string doesn't start/end with quotes and isn't a number, add quotes for MATLAB.
            // This assumes your MATLAB script expects these as string literals.
            if (typeof val === 'string' && !val.startsWith("'") && !val.endsWith("'") && isNaN(Number(val))) {
                formattedVal = `'${val}'`;
            }

            // Check if the variable assignment exists in the current code.
            // If it does, replace it. If not, the variable won't be "injected" into the displayed code string.
            // For truly dynamic injection of new variables, the backend's approach is superior.
            if (newCode.match(regex)) {
                newCode = newCode.replace(regex, `$1${key} = ${formattedVal}$2`);
            }
            // If you want to ADD variables that are NOT in the original file,
            // you would need more complex logic here (e.g., prepend them).
            // For now, we rely on the backend to inject ALL variables from the `variables` state.
        });
        setFileContent(newCode); // Update the content displayed in the <pre> tag
    }, [originalCode]); // Dependency: re-create if originalCode changes

    const runMatlabCode = () => {
    setIsLoading(true); // Start loading
    setMatlabResults(null); // Clear previous results

    // Ensure fileContent is updated with the very latest variables before sending.
    // This direct call ensures it's not debounced for the "Run" action.
    regenerateCodeWithVariables(variables); // Keep this line to update the displayed code if needed

    if (socketRef.current) {
        // ** CORRECTED: Pass the 'variables' object inside the 'args' array **
        socketRef.current.emit('run_matlab_code', {
            filename: fileName,
            // You can remove 'code: fileContent' if server.py runs the .m file by name and not by arbitrary code string.
            // If server.py has logic to execute arbitrary code string, then keep it.
            // For running existing file by name, 'code' field is not strictly necessary.
            args: [variables] // <--- THIS IS THE CRUCIAL CHANGE!
        });
        console.log("Emitting 'run_matlab_code' with payload:", {
            filename: fileName,
            // code: fileContent, // Optional: You can keep or remove this line based on server.py's implementation
            args: [variables] // Log the correct payload
        });
    } else {
        setMatlabResults({ error: "Socket not connected. Cannot run MATLAB code." });
        console.error("Socket.IO not connected!");
        setIsLoading(false); // Stop loading on error
    }
};

    // WebSocket setup for Flask backend and listeners
    useEffect(() => {
        socketRef.current = io('http://localhost:5000'); // Adjust URL if needed

        const onConnect = () => {
            setConnectionStatus("Connected");
            console.log('Connected to backend');
        };

        const onDisconnect = () => {
            setConnectionStatus("Disconnected");
            console.log('Disconnected from backend');
        };

        const onConnectError = (err) => {
            setConnectionStatus("Error");
            console.error('Connection error:', err);
        };

        // ** IMPORTANT: New Listener for 'matlab_code_output' from the backend **
        // This event carries the actual stdout/stderr from MATLAB execution.
        const handleMatlabCodeOutput = (data) => {
            setIsLoading(false); // Stop loading when output is received

            if (data.error) {
                setMatlabResults({ error: data.error, output: [] });
                setRunStatus("error");
            } else if (data.output) {
                setMatlabResults({ output: data.output, error: null });
                setRunStatus("success");
            } else {
                setMatlabResults({ output: ["No output or error received from MATLAB."], error: null });
                setRunStatus("success");
            }

        // Clear status after 3 seconds
        setTimeout(() => setRunStatus(""), 3000);
        };

        //Listener for matlab_result
        const handleMatlabResult = (data) => {
        setIsLoading(false);
        console.log("Received matlab_result:", data);

        if (data.error) {
            setMatlabResults({ error: data.error, output: [] });
            setRunStatus("error");
            toast.error("Error running MATLAB code.");
        } else {
            const combinedOutput = [
                ...(data.stdout || []),
                ...(data.results ? JSON.stringify(data.results, null, 2).split('\n') : [])
            ];
            setMatlabResults({ output: combinedOutput, error: null });
            setRunStatus("success");
            toast.success("MATLAB code finished running!");
        }

        setTimeout(() => setRunStatus(""), 3000);
    };

        // Register all listeners
        socketRef.current.on('connect', onConnect);
        socketRef.current.on('disconnect', onDisconnect);
        socketRef.current.on('connect_error', onConnectError);
        socketRef.current.on('matlab_code_output', handleMatlabCodeOutput); // Register the new listener
        socketRef.current.on('matlab_result', handleMatlabResult);


        // Cleanup function for useEffect
        return () => {
            if (socketRef.current) {
                socketRef.current.off('connect', onConnect);
                socketRef.current.off('disconnect', onDisconnect);
                socketRef.current.off('connect_error', onConnectError);
                socketRef.current.off('matlab_code_output', handleMatlabCodeOutput); // Unregister the new listener
                socketRef.current.off('matlab_result', handleMatlabResult);
                socketRef.current.disconnect();
            }
        };
    }, [setMatlabResults]); // Dependency: setMatlabResults prop

    // Initial regeneration when component mounts or originalCode/variables change
    // This ensures that when a file is uploaded or variables are changed, the displayed
    // code (`fileContent`) is immediately updated with the injected values.

    const downloadFile = () => {
        const blob = new Blob([fileContent], { type: "text/plain" });
        const link = document.createElement("a");
        link.href = URL.createObjectURL(blob);
        link.download = fileName.endsWith(".m") ? fileName : fileName + ".m";
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    };

    return (
        <div className={lightMode ? "light-mode" : ""}>
            <div className="main">
                <div className="content code-editor-layout">
                    <p>Connection Status: {connectionStatus}</p>
                    {/* Loading indicator for MATLAB execution */}
                    {isLoading && <p style={{ color: 'blue' }}>Running MATLAB code...</p>}

                    <div className="upload-section">
                        <h2>Upload File</h2>
                        <div style={{ display: "flex", gap: "10px", alignItems: "center" }}>
                            <label className="upload-btn">
                                <input
                                    type="file"
                                    id="matlabFile"
                                    hidden
                                    onChange={handleFileUpload}
                                    accept=".m,.txt"
                                />
                                <span className="icon">📂</span> Upload File
                            </label>
                            <button onClick={downloadFile} className="download-btn" disabled={!fileContent}>
                                ⬇️ Download
                            </button>
                        </div>
                        <p id="fileNameDisplay" className="file-name-display">
                            {fileName}
                        </p>
                    </div>

                    <div className="variables">
                        <h3>Variables</h3>
                        <div className="form-group">
                            {Object.entries(variables).map(([key, val]) => (
                                <React.Fragment key={key}>
                                    <label htmlFor={key}>{key}</label>
                                    <input
                                        id={key}
                                        type="text" // Keeping as text to allow flexible input for MATLAB strings/numbers
                                        value={val}
                                        onChange={(e) => {
                                            const newVal = e.target.value;
                                            setVariables((prev) => {
                                                const updated = { ...prev, [key]: newVal };
                                                regenerateCodeWithVariablesDebounced(updated); // Debounced update
                                                return updated;
                                            });
                                        }}
                                    />
                               </React.Fragment>
                            ))}
                        </div>
                    </div>

                    <div className="code-box">
                        <h2>Uploaded Code (with injected variables)</h2>
                        <textarea
                            id="codeEditor"
                            className="code-editor"
                            value={fileContent}
                            onChange={(e) => {
                                const updated = e.target.value;
                                setFileContent(updated);
                                localStorage.setItem("uploadedFileContent", updated);
                                setSaveStatus("Saved");

                                // Clear "Saved" after 1.5s
                                setTimeout(() => setSaveStatus(""), 1500);
                            }}
                            placeholder="// Uploaded code will appear here"
                            style={{
                                width: "100%",
                                height: "100%",
                                fontFamily: "monospace",
                                padding: "10px",
                                resize: "vertical",
                                backgroundColor: lightMode ? "#fff" : "#1e1e1e",
                                color: lightMode ? "#000" : "#eee",
                                border: "1px solid #ccc",
                                borderRadius: "5px",
                            }}
                        />
                        <button
                            className="save-button"
                            onClick={runMatlabCode}
                            disabled={!fileContent || isLoading || connectionStatus !== 'Connected'}
                        >
                            {isLoading ? 'Running...' : 'Run MATLAB Code'}
                        </button>

                        {runStatus === "success" && (
                            <p style={{ color: "green", fontWeight: "500", marginTop: "10px" }}>
                                ✅ MATLAB code finished running
                            </p>
                        )}
                        {runStatus === "error" && (
                        <p style={{ color: "red", fontWeight: "500", marginTop: "10px" }}>
                            ❌ Error running MATLAB code
                        </p>
                        )}

                        {saveStatus && (
                        <p style={{ color: "green", fontSize: "0.85rem", marginTop: "5px" }}>
                            {saveStatus}
                        </p>
                        )}
                    </div>
                </div>
            </div>
        </div>  
    );
}

export default CodeEditor;