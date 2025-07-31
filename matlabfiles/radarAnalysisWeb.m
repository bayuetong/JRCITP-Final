% filename: radarAnalysisWeb.m

function [graphs, results, reconstructedImageBase64] = radarAnalysisWeb(config, inputPaths)
% radarAnalysisWeb performs radar signal processing, breath/heartbeat detection,
% and image reconstruction, designed for integration with a web application.
%
% Inputs:
%   config: A struct containing analysis configuration parameters:
%     .runTests (boolean)
%     .sensor_test (boolean)
%     .hearbeat_test (boolean)
%     .speed_test (boolean)
%     .ts (scalar) - Sampling interval for JRC signals
%     .N_fft (scalar) - FFT size for PSD calculations
%     .desLen (scalar) - Desired length for frequency axis (usually N_fft)
%     .INUM (scalar) - Number of iterations for main loop (e.g., 500)
%     % IMPORTANT: Add any other scalar parameters your original script relies on here.
%     % For example, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U,
%     % M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, NP if they are
%     % meant to be configurable or are not hardcoded constants in your original script.
%
%   inputPaths: A struct with paths to required files, data, and scripts:
%     .mainScriptDir (char) - Directory where this MATLAB function (`radarAnalysisWeb.m`) resides.
%                             Also where any other `.m` helper functions (like `read_complex_binary.m`) are.
%     .usrpCaptureBaseDir (char) - Base directory for 'usrp_replay_captureX.dat' files.
%     .sensorCOM (char) - COM port for sensor (e.g., '/dev/ttyACM1'), if sensor_test is true.
%     .img_source_filename (char) - Original image filename (e.g., 'jarcat.jpg').
%                                   This should be a full path to the image on the server.
%     .rfnocTestScriptsDir (char) - Directory containing Python rfnoc_test scripts (e.g., image_rfnoc_replay.py).
%     .matlabImagesDir (char) - Directory containing Python img_encoder script (e.g., img_encoder.py).
%     .filterDataDir (char) - Directory for filter impulse response files (e.g., hf.dat, hf_hb.dat, hf_sensor.dat).
%
% Outputs:
%   graphs: A cell array of structs, where each struct represents a plot's data.
%     Each struct will have properties like:
%       .name (char) - Display name for the plot (matches dashboard.jsx's mapping)
%       .x_data (double array) - X-axis data
%       .y_data (double array) - Y-axis data
%       .z_data (double array, for heatmap like RDM) - Z-axis data for 3D plots
%       .type (char) - Plotly plot type (e.g., 'scatter' for line/points, 'heatmap')
%       .mode (char) - Plotly mode (e.g., 'lines', 'markers', 'lines+markers') for scatter plots
%       .y_type (char) - 'linear' or 'log' for y-axis scale (hint for frontend)
%       % You can add other Plotly-specific properties if needed, like 'line', 'marker', etc.
%
%   results: A struct with scalar analysis results (e.g., final breath rate, range, speed, BER).
%
%   reconstructedImageBase64: Base64 encoded string of the reconstructed image.
%     This allows embedding the image directly in HTML/JSON.

% --- 0. Initialize Output Variables ---
% It's good practice to pre-allocate or initialize all output variables.
graphs = {};
results = struct();
reconstructedImageBase64 = '';

% --- 1. Add Necessary Paths ---
% This is CRUCIAL. It tells MATLAB where to find your other `.m` files
% (like helper functions, or if your main logic is split) and where Python scripts are.
% Use `genpath` if your helper functions are in subdirectories within `mainScriptDir`.
try
    addpath(inputPaths.mainScriptDir);
    addpath(inputPaths.rfnocTestScriptsDir);
    addpath(inputPaths.matlabImagesDir);
    addpath(inputPaths.filterDataDir);
    fprintf('Paths added successfully.\n');
catch ME
    warning('Error adding paths: %s. Please ensure inputPaths fields are correct.', ME.message);
    results.error = ['Error adding MATLAB paths: ', ME.message];
    return; % Exit function early if paths cannot be added
end

% --- 2. Extract Configuration Parameters from `config` struct ---
% Assign input struct fields to local variables for easier use, mirroring your original script.
% Ensure all parameters your original script needs are listed here and
% provided by the Node.js `config` object.
runTests = config.runTests;
sensor_test = config.sensor_test;
hearbeat_test = config.hearbeat_test;
speed_test = config.speed_test;
ts = config.ts;
N_fft = config.N_fft;
desLen = config.desLen;
INUM = config.INUM; % Assuming this defines the length of your main processing loop

% --- 3. Extract File/Directory Paths from `inputPaths` struct ---
% These will be used to construct full file paths within your code.
usrpCaptureBaseDir = inputPaths.usrpCaptureBaseDir;
strCOM = inputPaths.sensorCOM;
img_source_filename = inputPaths.img_source_filename;


% --- 4. Suppress MATLAB Figure Windows (CRITICAL for Server Deployment) ---
% This prevents MATLAB from trying to open graphical windows on the server,
% which would cause errors or consume significant resources in a headless environment.
set(0, 'DefaultFigureVisible', 'off');
fprintf('MATLAB figure visibility set to OFF.\n');

% --- 5. Define Constants and Load Initial Data (from your original script) ---
% This section should contain all the initial setup from your original script.
% Constants, initial variable definitions, and loading of any necessary
% data files (like initial USRP captures, filter coefficients, etc.).

% Example Constants (Replace with actual values from your script):
cj = 1j; % Complex imaginary unit (built-in, but often defined for clarity)
ND = 1000; % Placeholder: Adjust to your actual value
NK = 100;  % Placeholder: Adjust to your actual value
MK = 20;   % Placeholder: Adjust to your actual value
NNP = 50;  % Placeholder: Adjust to your actual value
fc = 24e9; % Placeholder: Adjust to your actual value (e.g., carrier frequency)
N_antenna = 4; % Placeholder: Number of antennas
M = 256; % Placeholder: Number of samples in a frame
s = 1; % Placeholder: Scaling factor
B_data = 100e6; % Placeholder: Bandwidth
CP_len = 100; % Placeholder: Cyclic Prefix length
N_sqrt = 32; % Placeholder: Square root of samples
N = N_sqrt * N_sqrt; % Placeholder: Total samples
U = 16; % Placeholder: Upsampling factor
M_start = 1; % Placeholder: Starting index
M_max = 20; % Placeholder: Max M index
N_start = 1; % Placeholder: Starting N index
N_max = 20; % Placeholder: Max N index
MKNNP_samp_fc = 1e6; % Placeholder: Sampling frequency related constant
K1 = 1; % Placeholder: Another constant
samp_period = 1/MKNNP_samp_fc; % Placeholder: Sample period
NP = 1; % Placeholder: Another constant (e.g., number of frames)

% Load Filter Impulse Responses (ensure these .dat files are in `inputPaths.filterDataDir`)
try
    hf = load(fullfile(inputPaths.filterDataDir, 'hf_vs_RC_03_024_513_16_666HZ.dat'));
    hf_hb = load(fullfile(inputPaths.filterDataDir, 'hf_vs_RC_13_07_257_16_666HZ.dat'));
    hf_sensor = load(fullfile(inputPaths.filterDataDir, 'hf_vs_RC_sensor.dat'));
    fprintf('Filter impulse responses loaded.\n');
catch ME
    warning('Could not load filter files: %s. Check inputPaths.filterDataDir.', ME.message);
    results.error = ['Error loading filter files: ', ME.message];
    set(0, 'DefaultFigureVisible', 'on'); % Re-enable for debugging if exiting
    return;
end

% Initialize variables that will store collected data for plots (history)
% These will be populated inside your main processing loop.
angle_result = zeros(1, INUM); % Store raw angle data over time
speed_result = zeros(1, INUM); % Store speed estimates over time
range_result = zeros(1, INUM); % Store range estimates over time
respiratorTemp = []; % Will store raw sensor data if sensor_test is true

% Initial waveform setup (e.g., carrier, WaveImg1, pdsch, etc.)
% This is the part of your original script that sets up initial communication parameters
% or loads the first USRP capture for correlation/initial RDM.
% Replace these placeholders with your actual code.
% Example:
try
    % This assumes you have a `read_complex_binary` helper function
    v_initial_capture = read_complex_binary(fullfile(usrpCaptureBaseDir, 'usrp_replay_capture0.dat'));
    % Assuming WaveImg1 needs to be defined/loaded before resampling
    % You might need to load or generate 'WaveImg1' here if it's used before the img_propagation call.
    % Placeholder: WaveImg1 = randn(1, 1000); % Replace with actual WaveImg1 source
    % If WaveImg1 comes from img_propagation, consider how to initialize for the first signalAnalysis call.
    Orig1 = randn(1, 1228800); % Placeholder: Replace with actual source for Orig1
    % testMod = resample(v_initial_capture, 1228800, 2000000); % Example from your snippet
    % % ... Your correlation logic for `index` ...
    index = 100; % Placeholder for correlation index
    fprintf('Initial data loaded and variables initialized.\n');
catch ME
    warning('Error during initial data loading/setup: %s. Check usrpCaptureBaseDir and data files.', ME.message);
    results.error = ['Error in initial setup: ', ME.message];
    set(0, 'DefaultFigureVisible', 'on');
    return;
end


%% --- YOUR ORIGINAL RADAR ANALYSIS SCRIPT CONTENT GOES HERE ---
% IMPORTANT: Copy and paste the main logic of your original script here.
% AS YOU PASTE, MAKE THE FOLLOWING CRITICAL MODIFICATIONS:
%
% 1. REMOVE ALL PLOTTING COMMANDS:
%    - Delete `figure`, `plot`, `semilogy`, `imshow`, `nexttile`, `drawnow`, `animatedline`.
%    - Instead of plotting, ENSURE THE DATA that *would have been plotted* is stored
%      in variables (e.g., `x_data_for_plot1`, `y_data_for_plot1`, `z_data_for_rdm`, etc.).
%      These variables will be used later in Step 6 to populate the `graphs` output.
%
% 2. REPLACE HARDCODED FILE PATHS:
%    - Change paths like `'~/Desktop/...'` or `'/media/ramdisk/...'` to use `fullfile()`
%      and the `inputPaths` struct for portability.
%      - Example: `read_complex_binary('/media/ramdisk/usrp_replay_capture0.dat')` becomes
%                 `read_complex_binary(fullfile(usrpCaptureBaseDir, 'usrp_replay_capture0.dat'))`
%      - Example: `system('python3 ~/Desktop/rfnoc_test/image_rfnoc_replay.py')` becomes
%                 `system(['python3 ' fullfile(inputPaths.rfnocTestScriptsDir, 'image_rfnoc_replay.py')])`
%
% 3. ADAPT HELPER FUNCTION CALLS:
%    - If your original script calls helper functions (like `RunTest`, `signalAnalysis`, etc.)
%      which are now defined as local functions or separate `.m` files, ensure their
%      signatures match what you've defined (e.g., `RunTest` now needs `strCOM` and `rfnocTestScriptsDir`).
%
% 4. COLLECT HISTORY FOR "ANIMATED" PLOTS:
%    - If your original script used `animatedline` within a loop, you *must* modify that loop
%      to collect all data points into an array (e.g., `jrcBreathRatesHistory(loop_idx) = current_rate;`)
%      so that the entire time-series can be sent to the frontend for a static line plot.
%
% 5. IDENTIFY FINAL SCALAR RESULTS:
%    - Note down the variable names that hold your final scalar results (e.g., `p_rate_breath_JRC`, `range`, `v_est`, `ber`).
%      These will be used in Step 7.
%
% 6. ENSURE IMAGE RECONSTRUCTION OUTPUTS ARE READY:
%    - Make sure your image reconstruction part results in a `byteArray` variable (or similar)
%      containing the raw image data, and that it saves a temporary image file, if needed.
%      Also, ensure `ber` (Bit Error Rate) is calculated if you want to include it in results.

% --- Placeholder for Your Main Tests ---
% (e.g., calling RunTest, image_encoder, img_propagation)
% if runTests
%     fprintf('Running tests...\n');
%     if sensor_test
%         [~, respiratorTemp] = RunTest(strCOM, sensor_test, inputPaths.rfnocTestScriptsDir);
%     else
%         [~, ~] = RunTest('', sensor_test, inputPaths.rfnocTestScriptsDir); % Call without COM if no sensor
%     end
%     % Example: Call image_encoder. Ensure `img_source_filename` is a full path on the server.
%     image_encoder(img_source_filename, inputPaths.matlabImagesDir);
%     % Example: Call img_propagation. Ensure all its inputs are defined.
%     % [carrier, WaveImg1_updated, ...] = img_propagation(img_padded_length, img_input, NSlots, carrier, simLocal, ..., inputPaths.rfnocTestScriptsDir);
%     fprintf('Tests completed.\n');
% end

% --- Placeholder for Your Main Signal Processing Loop (e.g., JRC & Sensor analysis) ---
% This loop processes data for each time step, populating `angle_result`, `speed_result`, etc.
% And calculating breath/heartbeat rates over time for history plots.
fprintf('Starting main processing loop (%d iterations)...\n', INUM);
jrcBreathRatesHistory = zeros(1, INUM); % Example: To store rates over time
jrcHeartRatesHistory = zeros(1, INUM);
sensorBreathRatesHistory = zeros(1, INUM); % Only if sensor_test is true

for current_NUM = 1 : INUM
    % Example: Call signalAnalysis for each NUM
    % This helper function needs to be adapted to return all necessary data.
    % [current_angle, ~, ~, current_speed, current_range, x_ar_current, r_array_current, v_array_current] = ...
    %    signalAnalysis(current_NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP, usrpCaptureBaseDir);

    % Placeholder data for loop (REPLACE WITH YOUR ACTUAL CALCULATION)
    current_angle = sin(current_NUM * 0.1) + randn * 0.01;
    current_speed = abs(cos(current_NUM * 0.05));
    current_range = 5 + sin(current_NUM * 0.02);
    x_ar_current = randn(NNP, M_max); % Example RDM slice for current step
    r_array_current = linspace(0, 10, M_max);
    v_array_current = linspace(-1, 1, NNP);

    angle_result(current_NUM) = current_angle;
    speed_result(current_NUM) = current_speed;
    range_result(current_NUM) = current_range;

    % Example: Calculate current breath/heartbeat rates (REPLACE WITH YOUR ACTUAL CALCULATION)
    if current_NUM > 10 % Ensure enough data points for windowing
        % Apply your filtering, FFT, peak detection here using a window of `angle_result`
        % Example:
        % jrcBreathRatesHistory(current_NUM) = calculateJRCBreathRate(angle_result(max(1, current_NUM-len+1):current_NUM), ts, N_fft, desLen, hf);
        jrcBreathRatesHistory(current_NUM) = 10 + 2*sin(current_NUM * 0.03); % Dummy rate
        if hearbeat_test
            jrcHeartRatesHistory(current_NUM) = 60 + 5*cos(current_NUM * 0.02); % Dummy rate
        end
        if sensor_test && ~isempty(respiratorTemp) % Ensure sensor data exists
            % sensorBreathRatesHistory(current_NUM) = calculateSensorBreathRate(respiratorTemp(max(1, current_NUM_sensor_idx-len+1):current_NUM_sensor_idx), fs_sensor, N_fft, desLen, hf_sensor);
            sensorBreathRatesHistory(current_NUM) = 11 + sin(current_NUM * 0.04); % Dummy rate
        end
    end
end
fprintf('Main processing loop finished.\n');

% --- Placeholder for Your Image Reconstruction Logic ---
% This section processes data for image reconstruction and calculates BER.
% It should produce the `byteArray` of the reconstructed image and `ber`.
%
% Example:
byteArray = uint8(peaks(100) * 128 + 128); % Dummy byte array representing an image
ber = 0.01 + rand() * 0.005; % Dummy BER value
% Make sure your actual image data is in `byteArray` and `ber` is calculated.


% --- 6. Prepare `graphs` Output Data Structure for Frontend ---
% Transform the collected MATLAB data into the structured format expected by your frontend.
% The `name` property is crucial for `dashboard.jsx` to map data to the correct Plotly div.

% 6.1. JRC RDM (Range-Doppler Map) - Example uses the last calculated RDM or a combined one
% Assuming `x_ar_current`, `r_array_current`, `v_array_current` are available from the loop.
% Or if RDM is calculated only once at the beginning, use those variables.
graphs{end+1} = struct(...
    'name', 'JRC RDM', ...
    'rdm_data', struct('z', abs(x_ar_current').^2, 'x_axis', v_array_current, 'y_axis', r_array_current), ...
    'type', 'heatmap', ... % Explicitly set Plotly type
    'z_type', 'log' ... % Hint for frontend to use log scale for intensity
);

% 6.2. JRC Raw Breath Signal
graphs{end+1} = struct('name', 'JRC Raw Breath Signal', ...
                       'x_data', (ts:ts:INUM*ts)', ... % Time axis based on INUM and ts
                       'y_data', unwrap(angle_result)' - mean(unwrap(angle_result)), ... % Example data
                       'type', 'scatter', 'mode', 'lines', 'y_type', 'linear');

% 6.3. JRC Filtered Breath Signal (or Speed Result if speed_test is true)
if speed_test
    graphs{end+1} = struct('name', 'JRC Speed Result', ...
                           'x_data', (ts:ts:INUM*ts)', ...
                           'y_data', speed_result', ...
                           'type', 'scatter', 'mode', 'lines', 'y_type', 'linear');
else
    % You might need to re-filter the full `angle_result` here if `x2` isn't stored
    % from your original script's main loop.
    full_x1_unwrapped = unwrap(angle_result);
    full_x1_unwrapped = full_x1_unwrapped - mean(full_x1_unwrapped);
    x2_full_filtered = conv(full_x1_unwrapped, hf, 'same');
    graphs{end+1} = struct('name', 'JRC Filtered Breath Signal', ...
                           'x_data', (ts:ts:INUM*ts)', ...
                           'y_data', x2_full_filtered', ...
                           'type', 'scatter', 'mode', 'lines', 'y_type', 'linear');
end

% 6.4. JRC Breath PSD
% Assuming `spec_x1_full_breath_psd` and `x_axis_jrc_psd` are calculated from full signal
% You will need to calculate these from your `angle_result` or `x2_full_filtered` data.
spec_x1_full_breath_psd = fftshift(fft(x2_full_filtered, N_fft)) / sqrt(ND); % Example PSD calculation
x_axis_jrc_psd = (-desLen/2:desLen/2-1)*(1/ts)/desLen; % Example frequency axis
graphs{end+1} = struct('name', 'JRC Breath PSD', ...
                       'x_data', x_axis_jrc_psd', ...
                       'y_data', abs(spec_x1_full_breath_psd)'.^2, ...
                       'type', 'scatter', 'mode', 'lines', 'y_type', 'log');

% 6.5. JRC Heartbeat PSD (Conditional)
if hearbeat_test
    % Assuming `spec_x3_full_heart_psd` is calculated
    full_x3_filtered_heart = conv(full_x1_unwrapped, hf_hb, 'same'); % Use full unwrap signal
    spec_x3_full_heart_psd = fftshift(fft(full_x3_filtered_heart, N_fft)) / sqrt(ND);
    graphs{end+1} = struct('name', 'JRC Heartbeat PSD', ...
                           'x_data', x_axis_jrc_psd', ...
                           'y_data', abs(spec_x3_full_heart_psd)'.^2, ...
                           'type', 'scatter', 'mode', 'lines', 'y_type', 'log');
end

% 6.6. Sensor Breath PSD (Conditional)
if sensor_test && ~isempty(respiratorTemp)
    % Assuming `spec_x1_sensor_full_psd` is calculated from `respiratorTemp`
    spec_x1_sensor_full_psd = fftshift(fft(conv(respiratorTemp, hf_sensor, 'same'), N_fft)) / sqrt(ND);
    fs_sensor_psd = 10; % Example sensor sampling rate
    x_axis_sensor_psd = (-desLen/2:desLen/2-1)*fs_sensor_psd/desLen;
    graphs{end+1} = struct('name', 'Sensor Breath PSD', ...
                           'x_data', x_axis_sensor_psd', ...
                           'y_data', abs(spec_x1_sensor_full_psd)'.^2, ...
                           'type', 'scatter', 'mode', 'lines', 'y_type', 'log');
end

% 6.7. Respirator Monitor (Raw Sensor Data) (Conditional)
if sensor_test && ~isempty(respiratorTemp)
    graphs{end+1} = struct('name', 'Respirator Monitor', ...
                           'x_data', (0.1:0.1:length(respiratorTemp)*0.1)', ... % Adjust x-axis based on actual length
                           'y_data', respiratorTemp', ...
                           'type', 'scatter', 'mode', 'lines', 'y_type', 'linear');
end

% 6.8. Detected Rates Over Time (Multi-line plot)
% This requires `jrcBreathRatesHistory`, `jrcHeartRatesHistory`, `sensorBreathRatesHistory`
% to be populated with data from your main processing loop.
rate_plot_traces = {};
timeAxisRates = (1:INUM)' * ts; % Time axis for collected rates (adjust if your history is shorter)
rate_plot_traces{end+1} = struct('name', 'JRC Breath Rate', 'x_data', timeAxisRates, 'y_data', jrcBreathRatesHistory', 'mode', 'lines');
if sensor_test && ~isempty(sensorBreathRatesHistory)
    % Adjust timeAxisRates or ensure sensorBreathRatesHistory length matches
    rate_plot_traces{end+1} = struct('name', 'Respirator Sensor Rate', 'x_data', timeAxisRates, 'y_data', sensorBreathRatesHistory', 'mode', 'lines');
end
if hearbeat_test && ~isempty(jrcHeartRatesHistory)
    rate_plot_traces{end+1} = struct('name', 'JRC Heartbeat Rate', 'x_data', timeAxisRates, 'y_data', jrcHeartRatesHistory', 'mode', 'lines');
end
% This custom type tells your frontend to render multiple lines on one plot.
graphs{end+1} = struct('name', 'Detected Rates Over Time', ...
                       'traces', rate_plot_traces, ...
                       'type', 'multi-line-custom', ...
                       'x_label', 'Time (s)', ...
                       'y_label', 'Rate Per Minute');


% --- 7. Prepare `results` Output Struct (Scalar Analysis Results) ---
% Populate this struct with final calculated scalar values from your analysis.
% Make sure these variables are defined in your original script's logic.
% Example values are placeholders.
results.jrcBreathRate_final = jrcBreathRatesHistory(end); % Last calculated rate
results.finalRange = range_result(end); % Last calculated range
results.finalSpeed = speed_result(end); % Last calculated speed
results.ber = ber; % From image reconstruction

if hearbeat_test
    results.jrcHeartRate_final = jrcHeartRatesHistory(end);
end
if sensor_test && ~isempty(sensorBreathRatesHistory)
    results.sensorBreathRate_final = sensorBreathRatesHistory(end);
end


% --- 8. Handle Reconstructed Image Output (Base64 Encoding) ---
% This section converts your reconstructed image (from `byteArray`) into a Base64 string.
% Make sure `byteArray` contains your actual reconstructed image data (e.g., `uint8` values for a PNG/JPG).
% Determine the image extension based on what your reconstruction outputs.
img_extension = '.png'; % Example: '.png', '.jpg', '.bmp'
temp_img_filepath = fullfile(tempdir, ['reconstructed_img_' datestr(now,'yyyymmddHHMMSSFFF') img_extension]);

try
    % Ensure `byteArray` holds your actual reconstructed image data
    % Your original code's image reconstruction part should generate this.
    % fid = fopen(temp_img_filepath, 'wb');
    % fwrite(fid, byteArray, 'uint8'); % Write the image data to a temporary file
    % fclose(fid);

    % --- DUMMY IMAGE GENERATION (REMOVE THIS WHEN YOU INTEGRATE YOUR REAL CODE) ---
    % This is a placeholder to ensure the Base64 output works.
    % Replace this with your actual image saving logic from `byteArray`.
    dummy_img = peaks(100); % Example: some dummy 2D data
    dummy_img_scaled = uint8(255 * mat2gray(dummy_img)); % Scale to 0-255
    imwrite(dummy_img_scaled, temp_img_filepath, 'png'); % Save as PNG
    % --- END DUMMY IMAGE GENERATION ---

    if exist(temp_img_filepath, 'file')
        % Read the image data back from the temporary file
        % For PNG/JPEG, imread returns an M-by-N-by-3 array for color, M-by-N for grayscale.
        % For base64 encoding, you want the raw byte stream of the file.
        fid_read = fopen(temp_img_filepath, 'rb');
        img_bytes = fread(fid_read, inf, '*uint8');
        fclose(fid_read);

        reconstructedImageBase64 = matlab.net.base64encode(img_bytes); % Base64 encode the raw bytes
        delete(temp_img_filepath); % Clean up the temporary file
        fprintf('Reconstructed image encoded to Base64 and temporary file deleted.\n');
    else
        warning('Reconstructed image file not found at %s. Returning empty base64 string.', temp_img_filepath);
        reconstructedImageBase64 = '';
    end
catch ME
    warning('Error during image encoding: %s', ME.message);
    results.error = [results.error, newline, 'Image encoding error: ', ME.message];
    reconstructedImageBase64 = ''; % Ensure empty on error
end

% --- Cleanup: Re-enable Figure Visibility (Good Practice) ---
% This ensures MATLAB's default behavior is restored after this function call.
set(0, 'DefaultFigureVisible', 'on');
fprintf('MATLAB figure visibility set back to ON.\n');

end % End of the main radarAnalysisWeb function


%% --- Local Helper Functions (Place these at the end of this .m file) ---
% You MUST place the code for your helper functions here, or ensure they are
% in separate `.m` files located in `inputPaths.mainScriptDir`.
% Each function needs its own `function ... end` block.
%
% IMPORTANT: Modify these helper functions to:
%   1. Accept necessary paths as arguments (e.g., `usrpCaptureBaseDir`, `rfnocTestScriptsDir`).
%   2. Remove any plotting commands they might contain.
%   3. Return data explicitly if that data is needed by `radarAnalysisWeb`.

% Example structure for a helper function:
function v = read_complex_binary(filePath)
    % This is a placeholder. Insert your actual 'read_complex_binary' code here.
    % Example dummy implementation:
    % fid = fopen(filePath, 'rb');
    % if fid == -1
    %     error('Could not open file: %s', filePath);
    % end
    % v = fread(fid, inf, 'float32=>single');
    % v = v(1:2:end) + 1j*v(2:2:end);
    % fclose(fid);
    fprintf('Dummy read_complex_binary called for: %s\n', filePath);
    v = randn(1000, 1) + 1j * randn(1000, 1); % Dummy data
end

function [M_val, N_val] = RunTest(strCOM_path, sensor_flag, rfnocScriptsDir)
    % This is a placeholder. Insert your actual 'RunTest' code here.
    % It should call your Python scripts using `system()` with `fullfile()`.
    % Example:
    if sensor_flag
        fprintf('Running dummy sensor test on COM: %s\n', strCOM_path);
        % This part would typically call `sensorReading`
        N_val = 20 + 5*randn(500,1); % Dummy sensor data
    else
        N_val = [];
    end
    fprintf('Running dummy rfnoc_test from: %s\n', rfnocScriptsDir);
    % system(['python3 ' fullfile(rfnocScriptsDir, 'rfnoc_replay_30k.py')]); % Actual call
    M_val = -1; % Dummy value
end

function tmp = sensorReading(strCOM)
    % This is a placeholder. Insert your actual 'sensorReading' code here.
    % Remove any `plot` or `figure` commands.
    fprintf('Dummy sensorReading called on COM: %s\n', strCOM);
    tmp = 25 + 5*sin(linspace(0, 10*pi, 100))' + randn(100,1); % Dummy temperature data
end

function [angle_res, angle_res1, angle_res2, speed_res, range_res, x_ar_out, r_arr_out, v_arr_out] = ...
         signalAnalysis(NUM, Orig1_in, index_in, NK_in, MK_in, NNP_in, fc_in, N_antenna_in, M_in, s_in, B_data_in, CP_len_in, N_sqrt_in, N_in, U_in, M_start_in, M_max_in, N_start_in, N_max_in, MKNNP_samp_fc_in, K1_in, samp_period_in, cj_in, NP_in, usrpCaptureBaseDir_in)
    % This is a placeholder. Insert your actual 'signalAnalysis' code here.
    % Ensure it takes all necessary inputs and returns the specified outputs.
    % Replace file loads with `fullfile(usrpCaptureBaseDir_in, ...)`
    % Ensure `x_ar_out`, `r_arr_out`, `v_arr_out` for RDM are calculated and returned.
    fprintf('Dummy signalAnalysis called for NUM=%d\n', NUM);
    angle_res = rand() * 10;
    angle_res1 = rand() * 5;
    angle_res2 = rand() * 5;
    speed_res = rand() * 0.5;
    range_res = rand() * 10;

    % Dummy RDM data (replace with actual calculation)
    x_ar_out = randn(NNP_in, M_max_in);
    r_arr_out = linspace(0, 10, M_max_in);
    v_arr_out = linspace(-1, 1, NNP_in);
end

function image_encoder(img_filename_input, matlabImagesDir_in)
    % This is a placeholder. Insert your actual 'image_encoder' code here.
    % Replace system call: `system(['python3 ' fullfile(matlabImagesDir_in, 'img_encoder.py') ' ' img_filename_input]);`
    fprintf('Dummy image_encoder called for %s from %s\n', img_filename_input, matlabImagesDir_in);
end

function [carrier_out, WaveImg1_out, dmrsAntSymbols_out, dmrsAntIndices_out, dmrsSymbols_all_out, dmrsIndices_all_out, trBlk_tx_out, trBlkSizes_out, pdsch_all_out, pdschIndices_out, pdschIndices_all_out, encodeDLSCH_out] = ...
         img_propagation(img_padded_length_in, img_input_in, NSlots_in, carrier_in, simLocal_in, pdsch_in, encodeDLSCH_in, decodeDLSCHLocal_in, pdschextra_in, harqEntity_in, rfnocTestScriptsDir_in)
    % This is a placeholder. Insert your actual 'img_propagation' code here.
    % Replace save calls with `save(fullfile(tempdir, ...))` if temporary.
    % Replace system call: `system(['python3 ' fullfile(rfnocTestScriptsDir_in, 'image_rfnoc_replay.py')]);`
    fprintf('Dummy img_propagation called from %s\n', rfnocTestScriptsDir_in);

    % Dummy outputs
    carrier_out = struct('SubcarrierSpacing', 15); % Example dummy struct
    WaveImg1_out = randn(1, 100);
    dmrsAntSymbols_out = []; dmrsAntIndices_out = []; dmrsSymbols_all_out = []; dmrsIndices_all_out = [];
    trBlk_tx_out = []; trBlkSizes_out = []; pdsch_all_out = []; pdschIndices_out = []; pdschIndices_all_out = [];
    encodeDLSCH_out = [];
end