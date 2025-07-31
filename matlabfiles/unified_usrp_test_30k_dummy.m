% unified_usrp_test_30k_dummy.m
function [graphs, results] = unified_usrp_test_30k_dummy()
% This function generates dummy radar analysis data and prepares it for
% display in the frontend. It does not open MATLAB figures.

% Dummy Data Setup
INUM = 1000;
ts = 0.06;
angle_result = randn(1, INUM);        % Random phase data
% Ensure angle_result is not all zeros or too small, to avoid issues with unwrap
if all(angle_result == 0)
    angle_result = randn(1,INUM) + 0.1;
end
speed_result = rand(1, INUM);         % Random speed
respiratorTemp = rand(1, INUM);       % Random respiration sensor data
hf = ones(1, 10) / 10;                % Simple moving average filter
hf_hb = ones(1, 5) / 5;               % Heartbeat filter
hf_sensor = ones(1, 10) / 10;         % Sensor filter

range_dummy = rand * 5;               % Dummy range in meters
v_est_dummy = rand * 2;               % Dummy velocity estimate
hearbeat_test = true;                 % Keep these as internal flags for dummy data
sensor_test = true;                   % Keep these as internal flags for dummy data
speed_test = true;                    % Keep these as internal flags for dummy data


% Process signals
x1 = unwrap(angle_result);
x1 = x1 - mean(x1);
x2 = conv(x1, hf, 'same');
x3 = conv(x1, hf_hb, 'same');

% Time axis for plotting
time_axis = (ts:ts:INUM*ts);

% FFT and PSD for Breath Signal
N_fft = 8192;
fs = 1 / ts;
spec_x1 = fftshift(fft(x2, N_fft)) / sqrt(INUM);
f_axis = (-N_fft/2:N_fft/2-1) * fs / N_fft;

% Estimate breath rate (JRC)
f_start = 0.1;
MS = floor(f_start / fs * N_fft);
y_pf = abs(spec_x1(N_fft/2 + MS:N_fft)).^2;
[~, locs] = max(y_pf);
p_freq_breath_JRC = (locs + MS - 2) * fs / N_fft;
p_rate_breath_JRC = p_freq_breath_JRC * 60;

% Optional: Heartbeat estimation (JRC)
p_rate_heart_JRC = NaN; % Initialize
if hearbeat_test
    spec_x3 = fftshift(fft(x3, N_fft)) / sqrt(INUM);
    y_pf_hb = abs(spec_x3(N_fft/2 + MS:N_fft)).^2;
    [~, locs_hb] = max(y_pf_hb);
    p_freq_heart_JRC = (locs_hb + MS - 2) * fs / N_fft;
    p_rate_heart_JRC = p_freq_heart_JRC * 60;
end

% Optional: Sensor data processing
p_rate_breath_Sensor = NaN; % Initialize
if sensor_test
    x2_sensor = conv(respiratorTemp, hf_sensor, 'same');
    spec_x_sensor = fftshift(fft(x2_sensor, N_fft)) / sqrt(INUM);
    fs_sensor = 10; % Assuming a sample rate for the sensor data
    y_pf_sensor = abs(spec_x_sensor(N_fft/2 + MS:N_fft)).^2;
    [~, locs_s] = max(y_pf_sensor);
    p_freq_sensor = (locs_s + MS - 2) * fs_sensor / N_fft;
    p_rate_breath_Sensor = p_freq_sensor * 60;
end

% Prepare graphs output for frontend
graphs = {}; % Initialize as a cell array

% Raw Breath Signal
graphs{end+1} = struct(...
    'name', 'JRC Raw Breath', ...
    'x_data', time_axis', ... % Transpose to column vector
    'y_data', x1' ...         % Transpose to column vector
);

% Filtered Breath Signal
graphs{end+1} = struct(...
    'name', 'JRC Filtered Breath', ...
    'x_data', time_axis', ... % Transpose to column vector
    'y_data', x2' ...         % Transpose to column vector
);

% Breath Signal Spectrum (PSD)
graphs{end+1} = struct(...
    'name', 'JRC Breath PSD', ...
    'x_data', f_axis', ...        % Transpose to column vector
    'y_data', abs(spec_x1).^2' ... % Transpose to column vector
);

% Heartbeat PSD (if enabled)
if hearbeat_test && ~isnan(p_rate_heart_JRC)
    graphs{end+1} = struct(...
        'name', 'JRC Heartbeat PSD', ...
        'x_data', f_axis', ...        % Transpose to column vector
        'y_data', abs(spec_x3).^2' ... % Transpose to column vector
    );
end

% Prepare scalar results output for frontend
results = struct(...
    'jrcBreathRate', p_rate_breath_JRC, ...
    'jrcHeartRate', p_rate_heart_JRC, ...
    'sensorBreathRate', p_rate_breath_Sensor, ...
    'range', range_dummy, ...
    'speed', v_est_dummy ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User-defined Functions' Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These functions are kept as local functions within this file.

function [M, N] = RunTest(strCOM,sensor_test)
% Dummy RunTest function
M = rand(10,1); % Dummy output
N = rand(10,1); % Dummy output
disp("Dummy Tests completed!");

function tmp = sensorReading(strCOM)
% Dummy sensorReading function
tmp = rand(600,1); % Dummy output
disp("Dummy sensor data collected.");

function [angle_result,angle_result_1,angle_result_2,speed_result,range_result]= signalAnalysis(NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP)
% Dummy signalAnalysis function
angle_result = rand;
angle_result_1 = rand;
angle_result_2 = rand;
speed_result = rand;
range_result = rand;
disp("Dummy signal analysis completed.");

function [carrier, WaveImg1, dmrsAntSymbols, dmrsAntIndices, dmrsSymbols_all, dmrsIndices_all, trBlk_tx, trBlkSizes, pdsch_all, pdschIndices, pdschIndices_all, encodeDLSCH] = img_propagation(img_padded_length, img_input, NSlots, carrier, simLocal, pdsch, encodeDLSCH, decodeDLSCHLocal, pdschextra, harqEntity)
% Dummy img_propagation function
WaveImg1 = rand(100,1); % Dummy output
carrier = struct('NSlot', 0); % Dummy carrier
dmrsAntSymbols = rand(10,1); dmrsAntIndices = randi(100,10,1);
dmrsSymbols_all = rand(10,1); dmrsIndices_all = randi(100,10,1);
trBlk_tx = randi([0,1],10,1); trBlkSizes = 10;
pdsch_all = {}; pdschIndices = randi(100,10,1); pdschIndices_all = randi(100,10,1);
encodeDLSCH = []; % Dummy handle

disp('Dummy image propagation...');
% Remove system calls if not needed for dummy test
% system('python3 ~/Desktop/rfnoc_test/image_rfnoc_replay.py');

function a = rfnoc_test(~)
% Dummy rfnoc_test function
disp("Dummy rfnoc_test.py called.");
% Remove system calls if not needed for dummy test
% system('python3 ~/Desktop/rfnoc_test/rfnoc_replay_30k.py');
a = -1;

function image_encoder(img_filename)
% Dummy image_encoder function
disp('Dummy image encoding...');
% Remove system calls if not needed for dummy test
% system(['python3 ~/Desktop/matlab/images/img_encoder.py ' img_filename]);

