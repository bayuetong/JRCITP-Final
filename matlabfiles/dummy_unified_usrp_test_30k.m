runTests=1;
sensor_test=true;
hearbeat_test=true;
speed_test=false;
real_time=true;
ts=0.06;

% IMPORTANT: Define function signature to accept arguments from Python
function unified_usrp_test_30k(runTests, sensor_test, hearbeat_test, speed_test, real_time, ts, img_filename)

p = gcp('nocreate');
if (isempty(p)); parpool(6); end % Adjusted parpool size to 6 based on your provided file

addpath("~/Desktop/matlab/Sensor_code"); addpath("~/Desktop/matlab/test_data");
addpath('/home/xzhang/Documents/MATLAB/Examples/R2023b/5g/NewRadioPUSCHThroughputExample'); addpath('~/Desktop/matlab/images');
addAttachedFiles(gcp,["read_complex_binary.m" "ProcessBuffer.m"])

%% Tests Variables %%
strCOM = ['/dev/ttyACM1' ...
    ''];

% Dummy load for filter coefficients
hf=rand(1, 100);
hf_hb=rand(1, 100);
hf_sensor=rand(1, 100);
hf_25HZ=rand(1, 100);
hf_33HZ=rand(1, 100);
hf_hb_25HZ=rand(1, 100);
hf_hb_33HZ=rand(1, 100);

if ts==0.03
    hf=hf_33HZ;
    hf_hb=hf_hb_33HZ;
elseif ts==0.04
    hf=hf_25HZ;
    hf_hb=hf_hb_25HZ;
end

%% USRP signal analysis variables initialisation %%
% Dummy load for waveform, waveforminfo, and signal
waveformInfo.Nfft = 4096;
waveformInfo.CyclicPrefixLengths = ones(1,14) * 288; % Example CP lengths
waveformInfo.SampleRate = 200000; % Example sample rate (Hz)
waveformInfo.NumSubframes = 10; % Example number of subframes
waveformInfo.SlotsPerSubframe = 1; % Example slots per subframe
waveformInfo.SymbolsPerSlot = 14; % Example symbols per slot

% Dummy initialization of parameters needed by signalAnalysis_dummy
% (These would normally come from actual setup or previous calculations)
Orig1 = randn(100,1); % Dummy data
index = 1:100;
NK = 10; MK = 5; NNP = 2; fc = 2.4e9; % Example values
N_antenna = 4; M = 10; s = 1; B_data = rand(10,1); CP_len = 5;
N_sqrt = 2; N = 20; U = 1; M_start = 1; M_max = 10; N_start = 1; N_max = 20;
MKNNP_samp_fc = 100; K1 = 5; samp_period = 0.01; cj = 1; NP = 10;

% Initializing result arrays for collection after parfor
angle_result = zeros(1, 1000); % Assuming INUM up to 1000 for this example
angle_result_1 = zeros(1, 1000);
angle_result_2 = zeros(1, 1000);
range_result = zeros(1, 1000);
speed_result = zeros(1, 1000);

INUM = 1000; % Example: Loop 1000 times

parfor NUM = 1:INUM
    % Call your signal analysis function (using dummy for demonstration)
    % [angle,angle_1,angle_2,speed,range_val]=signalAnalysis_dummy(NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP);

    % --- Dummy real-time data generation for demonstration ---
    % Replace this with actual outputs from your signalAnalysis_dummy or equivalent logic
    current_angle = rand() * 2 * pi;
    current_angle_1 = randn();
    current_angle_2 = randn();
    current_speed = randn() * 10;
    current_range_val = randn() * 10;

    % Dummy data for frame-based plots
    time_points_in_frame = 1:100; % Example: 100 samples per frame
    current_breath_signal = sin(2*pi*0.5*time_points_in_frame/100) + randn(1,100)*0.1; % Example sine wave + noise
    current_breath_psd = abs(fft(current_breath_signal)); % Simple PSD for dummy
    current_breath_psd = current_breath_psd(1:floor(end/2)); % Take first half for meaningful PSD
    
    % Dummy RDM data (e.g., 20x20 matrix)
    rdm_rows = 20; rdm_cols = 20;
    current_rdm_data = rand(rdm_rows, rdm_cols);
    current_rdm_x_array = linspace(-10, 10, rdm_cols); % Dummy speed axis
    current_rdm_y_array = linspace(0, 5, rdm_rows);   % Dummy range axis

    % Store results (if needed for final collection after loop)
    angle_result(NUM) = current_angle;
    angle_result_1(NUM) = current_angle_1;
    angle_result_2(NUM) = current_angle_2;
    range_result(NUM) = current_range_val;
    speed_result(NUM) = current_speed;

    % Create a struct to send data for the current iteration
    % Ensure field names match what dashboard.jsx expects ('timestamp', 'angle', 'speed', 'range', etc.)
    data_to_send = struct(...
        'timestamp', now(), ... % Use MATLAB's current time (datenum) or simple counter for x-axis
        'angle', current_angle, ...
        'speed', current_speed, ...
        'range', current_range_val, ...
        'x1_segment', current_breath_signal, ...
        'spec_x1_segment', current_breath_psd, ...
        'rdm_data', current_rdm_data, ...
        'x_array', current_rdm_x_array, ...
        'y_array', current_rdm_y_array ...
    );

    % --- NEW: Call the Python function to emit real-time data ---
    % 'eng' is your MATLAB engine instance, which now has access to 'emit_realtime_data_python'
    % nargout=0 indicates that this call does not return any output to MATLAB.
    try
        eng.emit_realtime_data_python(data_to_send, nargout=0);
    catch ME
        warning('Failed to emit real-time data to Python: %s', ME.message);
    end
    
    % Use drawnow to ensure graphics/output are processed,
    % though for engine calls, it mainly applies to MATLAB's own figures.
    % drawnow;
    
    % Display message (will be out of sequence due to parfor, as discussed)
    disp("parfor" + string(NUM));

    % Optional: Add a small pause to simulate processing time and control update rate
    % pause(0.01); % Pause for 10 milliseconds
end

% You can emit final results here if your script processes them all at the end
final_results_struct = struct(...
    'JRC Raw Breath', struct('data', rand(1,100), 'time_axis', 1:100, 'type', 'line', 'name', 'JRC Raw Breath'), ...
    'JRC Filtered Breath', struct('data', rand(1,100), 'time_axis', 1:100, 'type', 'line', 'name', 'JRC Filtered Breath'), ...
    'JRC Breath PSD', struct('data', rand(1,50), 'freq_axis', linspace(0, 10, 50), 'type', 'line-logy', 'name', 'JRC Breath PSD'), ...
    'JRC Heartbeat PSD', struct('data', rand(1,50), 'freq_axis', linspace(0, 10, 50), 'type', 'line-logy', 'name', 'JRC Heartbeat PSD'), ...
    'JRC RDM', struct('data', rand(20,20), 'x_array', linspace(-10,10,20), 'y_array', linspace(0,5,20), 'type', 'mesh3d', 'name', 'JRC RDM') ...
);

% Assuming a function in Python to emit final data as well
eng.emit_final_data_python(final_results_struct, nargout=0);
eng.simulation_complete_python(struct('message', 'MATLAB simulation finished successfully!'), nargout=0);

disp("MATLAB script finished execution.");

end % End of function

% (You might need signalAnalysis_dummy.m in your matlabfiles directory or defined within this script if not already)
% function [angle,angle_1,angle_2,speed,range_val]=signalAnalysis_dummy(NUM, Orig1, index, NK, MK, NNP, fc, N_antenna, M, s, B_data, CP_len, N_sqrt, N, U, M_start, M_max, N_start, N_max, MKNNP_samp_fc, K1, samp_period, cj, NP)
%     % This is a dummy implementation based on your code snippet.
%     % Replace with your actual signal analysis logic.
%     angle = rand() * 2 * pi;
%     angle_1 = randn();
%     angle_2 = randn();
%     speed = randn() * 10;
%     range_val = randn() * 10;
% end
