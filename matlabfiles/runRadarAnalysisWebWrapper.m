% runRadarAnalysisWebWrapper.m

function [graphs, results, reconstructedImageBase64] = runRadarAnalysisWebWrapper(config, inputPaths)
% This wrapper script executes the original radar analysis code
% and captures its outputs for web visualization.

% --- 0. Initialize Outputs ---
graphs = {};
results = struct();
reconstructedImageBase64 = '';

% --- 1. Add Necessary Paths ---
% Ensure MATLAB can find your original script, helper functions,
% and the Python scripts if they're in different directories.
addpath(inputPaths.originalScriptDir); % Directory where original_radar_analysis_script.m is
addpath(inputPaths.matlabHelperScriptsDir); % e.g., read_complex_binary, filter data files
addpath(inputPaths.rfnocTestScriptsDir);   % Python scripts called by system()
addpath(inputPaths.matlabImagesDir);      % Python image encoder script

% --- 2. Set Up Global Variables / Simulate External Inputs ---
% Your original script might rely on variables being defined in the workspace
% (e.g., `ts`, `INUM`, `speed_test`, `respiratorTemp`, `Orig1`, etc.).
% These need to be passed from Node.js via `config` and `inputPaths` structs
% and then made available to the original script.

% Directly assign from config struct
ts = config.ts;
INUM = config.INUM; % Assuming INUM is part of config or derived
speed_test = config.speed_test;
hearbeat_test = config.hearbeat_test;
sensor_test = config.sensor_test;
N_fft = config.N_fft;
desLen = config.desLen;

% File paths for the original script (modify if your script expects different variable names)
usrpCaptureBaseDir = inputPaths.usrpCaptureBaseDir;
strCOM = inputPaths.sensorCOM; % Will be passed to sensorReading/RunTest
img_filename = inputPaths.img_source_filename; % For image encoding/decoding

% --- Simulate `load` commands from your original script ---
% If your original script `load`s .dat files, make sure they are in `matlabHelperScriptsDir`
% or update the `load` commands within your original script to use `fullfile(inputPaths.matlabHelperScriptsDir, 'filename.dat')`.
% This is one of the very few *minor* changes you might need in the original script.
% If not, you might have to load them *here* in the wrapper and make them global, which is less clean.
% Example (if original script's `load`s are not adjusted for paths):
% global hf hf_hb hf_sensor; % Make them global if original script accesses them as globals
% hf = load(fullfile(inputPaths.matlabHelperScriptsDir, 'hf_vs_RC_03_024_513_16_666HZ.dat'));
% hf_hb = load(fullfile(inputPaths.matlabHelperScriptsDir, 'hf_vs_RC_13_07_257_16_666HZ.dat'));
% hf_sensor = load(fullfile(inputPaths.matlabHelperScriptsDir, 'hf_vs_RC_sensor.dat'));

% --- IMPORTANT: Suppress figure windows in MATLAB engine ---
% This prevents MATLAB from trying to open GUI windows on the server.
% This is critical for headless server environments.
set(0, 'DefaultFigureVisible', 'off');

% --- 3. Execute the Original Script ---
% The `run` command executes the script in the current workspace.
% This is the magic step!
try
    run(fullfile(inputPaths.originalScriptDir, 'original_radar_analysis_script.m'));
catch ME
    % Re-enable figure visibility in case of error
    set(0, 'DefaultFigureVisible', 'on');
    rethrow(ME); % Rethrow the error to Node.js
end

% Re-enable figure visibility (important if you run other things later or for other calls)
set(0, 'DefaultFigureVisible', 'on');


% --- 4. Capture Data from Workspace Variables ---
% After `run`, all variables created by `original_radar_analysis_script.m`
% will be present in the workspace of `runRadarAnalysisWebWrapper.m`.
% Now, extract the data you want to send to the web.

% Example: Capture data for graphs
% You'll need to know the exact variable names your original script uses for plot data.
% For example, `x1`, `x2`, `spec_x1`, `respiratorTemp`, `angle_result`, `speed_result`, `x_ar`, `r_array`, `v_array`, etc.

% JRC Raw Breath Signal
graphs{end+1} = struct('name', 'JRC Raw Breath Signal', ...
                       'x_data', (ts:ts:INUM*ts)', ...
                       'y_data', x1', ... % Assuming x1 is available after running original script
                       'type', 'line', 'y_type', 'linear');

% JRC Speed Result / JRC breath signal after a filter
if speed_test
    graphs{end+1} = struct('name', 'JRC Speed Result', ...
                           'x_data', (ts:ts:INUM*ts)', ...
                           'y_data', speed_result', ... % Assuming speed_result is available
                           'type', 'line', 'y_type', 'linear');
else
    graphs{end+1} = struct('name', 'JRC Filtered Breath Signal', ...
                           'x_data', (ts:ts:INUM*ts)', ...
                           'y_data', x2', ... % Assuming x2 is available
                           'type', 'line', 'y_type', 'linear');
end

% JRC Breath Signal PSD
x_axis_jrc = (-desLen/2:desLen/2-1)*(1/ts)/desLen; % Reconstruct x_axis if not explicitly stored
graphs{end+1} = struct('name', 'JRC Breath Signal PSD', ...
                       'x_data', x_axis_jrc', ...
                       'y_data', abs(spec_x1)'.^2, ... % Assuming spec_x1 is available
                       'type', 'line', 'y_type', 'log');

% JRC Heartbeat Signal PSD
if hearbeat_test
    % Assuming spec_x3 is available from original script after hearbeat_test is true
    graphs{end+1} = struct('name', 'JRC Heartbeat Signal PSD', ...
                           'x_data', x_axis_jrc', ...
                           'y_data', abs(spec_x3)'.^2, ...
                           'type', 'line', 'y_type', 'log');
end

% Sensor Breath Signal PSD
if sensor_test
    fs_sensor_psd = 10; % Assuming this is fixed for sensor
    x_axis_sensor_psd = (-desLen/2:desLen/2-1)*fs_sensor_psd/desLen;
    % Assuming spec_x1_sensor is available
    graphs{end+1} = struct('name', 'Sensor Breath Signal PSD', ...
                           'x_data', x_axis_sensor_psd', ...
                           'y_data', abs(spec_x1_sensor)'.^2, ...
                           'type', 'line', 'y_type', 'log');
end

% Respirator Monitor (Raw Sensor Data)
if sensor_test
    % Assuming respiratorTemp is available and correctly sized
    graphs{end+1} = struct('name', 'Respiration Monitor (Raw Sensor)', ...
                           'x_data', (0.1:0.1:length(respiratorTemp)*0.1)', ...
                           'y_data', respiratorTemp', ...
                           'type', 'line', 'y_type', 'linear');
end

% JRC RDM (Range-Doppler Map)
% Assuming x_ar, r_array, v_array are available after running original script
graphs{end+1} = struct(...
    'name', 'JRC RDM', ...
    'x_data', v_array, ...
    'y_data', r_array(1:min(20, length(r_array))), ... % Only first 20 points
    'z_data', abs(x_ar(:, 1:min(20, length(r_array)))').^2, ...
    'type', 'heatmap', ...
    'z_type', 'log' ...
);


% Combined Rates Over Time (This is trickier as your original script uses animatedline)
% The wrapper cannot "see" the history of `addpoints` directly.
% **MAJOR CAVEAT:** If your original script *only* outputs `p_rate_breath_JRC` as a single scalar
% at the *end* of the script, then you can only get the final rate.
% If your original script calculates these rates within its loops and stores them in arrays (e.g., `jrcBreathRates`),
% *then* you can capture those arrays here. If not, the "animated" plot is very hard to reconstruct.
% You might have to modify your original script slightly to store these history arrays.
% Let's assume for now your original script *does* create final scalar variables like `p_rate_breath_JRC`.
% If you need the historical animated line, you *must* change the original code to collect these into arrays.

% Example for final scalar results (assuming original script produces these)
results.jrcBreathRate_final = p_rate_breath_JRC; % Assuming this variable exists
if hearbeat_test
    results.jrcHeartRate_final = p_rate_heart_JRC;
end
if sensor_test
    results.sensorBreathRate_final = p_rate_breath_Sensor;
end
results.finalRange = range; % Assuming 'range' variable holds the final value
results.finalSpeed = v_est; % Assuming 'v_est' variable holds the final value
results.ber = ber; % Assuming 'ber' is calculated in image section


% --- 5. Handle Reconstructed Image Output ---
% The original script saves the image to `~/Desktop/matlab/test_data/temp.ext` and then displays it.
% We need to capture that saved file.
img_extension = img_filename(regexp(img_filename, '\.'): end); % From original script's logic
temp_img_path = fullfile('~/Desktop/matlab/test_data/', ['temp' img_extension]); % This path is hardcoded in original script

% You need to ensure the original script has permissions to write to this path on the server.
% Better: modify your original script to save to `tempdir` using `fullfile(tempdir, ...)`
% If not, you must ensure the hardcoded path exists and is writable on the server.

if exist(temp_img_path, 'file')
    img_data = fileread(temp_img_path);
    reconstructedImageBase64 = matlab.net.base64encode(uint8(img_data));
    delete(temp_img_path); % Clean up the temporary file
else
    warning('Reconstructed image file not found at %s', temp_img_path);
    reconstructedImageBase64 = ''; % Return empty if not found
end

% --- Cleanup: Remove paths added by this wrapper if desired ---
% rmpath(inputPaths.originalScriptDir);
% rmpath(inputPaths.matlabHelperScriptsDir);
% rmpath(inputPaths.rfnocTestScriptsDir);
% rmpath(inputPaths.matlabImagesDir);

end % End of wrapper function