% filename: runAnalysisWrapper.m
% Path: runAnalysisWrapper.m

function [output_graph_data, output_analysis_results, errors] = runAnalysisWrapper(input_config_struct)
% This function now generates dummy graph and analysis data for frontend testing.
% It no longer relies on external MATLAB scripts or data files.
% The output structure is adjusted to match dashboard.jsx's createPlotlyTrace expectations.

% Initialize outputs
output_graph_data = {};
output_analysis_results = struct();
errors = '';

try
    % --- 1. Process Input Configuration ---
    % These inputs are acknowledged and can be used to vary dummy data if needed.
    runTests = input_config_struct.runTests;
    sensor_test = input_config_struct.sensor_test;
    hearbeat_test = input_config_struct.hearbeat_test;
    speed_test = input_config_struct.speed_test;
    real_time = input_config_struct.real_time;
    ts = input_config_struct.ts;
    img_filename = input_config_struct.img_filename;
    strCOM = input_config_struct.strCOM;

    fprintf('runAnalysisWrapper (DUMMY MODE) received config:\n');
    disp(input_config_struct);
    fprintf('Generating dummy data...\n');

    % --- 2. Generate Dummy Graph Data (matching dashboard.jsx's createPlotlyTrace structure) ---
    % Ensure 'name' property matches the graphNameMapping in dashboard.jsx

    % Common dummy data for line plots
    num_points = 200;
    % Ensure these are column vectors for better JSON serialization consistency
    dummy_time = (0:num_points-1)' * 0.05; % Use transpose ' to make it a column vector
    dummy_freq_psd = linspace(0, 5, num_points)'; % Use transpose ' to make it a column vector

    % 1. Dummy JRC RDM (Heatmap)
    % This is the 3D data, expecting 'rdm_data.z', 'rdm_data.x_axis', 'rdm_data.y_axis'
    num_rdm_rows = 50;
    num_rdm_cols = 40;
    dummy_z_rdm = rand(num_rdm_rows, num_rdm_cols) * 100; % A 2D matrix
    % Add a prominent peak to make it visually clear
    dummy_z_rdm(round(num_rdm_rows/3), round(num_rdm_cols/2)) = 150;
    dummy_z_rdm(round(num_rdm_rows/2), round(num_rdm_cols/4)) = 120;


    dummy_x_rdm_axis = linspace(-5, 5, num_rdm_cols); % Dummy Doppler axis (row vector is fine for axes)
    dummy_y_rdm_axis = linspace(0, 10, num_rdm_rows); % Dummy Range axis (row vector is fine for axes)

    rdm_data = struct();
    rdm_data.name = 'JRC RDM';
    rdm_data.rdm_data = struct(); % Nested struct as expected by JS
    rdm_data.rdm_data.z = dummy_z_rdm;
    rdm_data.rdm_data.x_axis = dummy_x_rdm_axis;
    rdm_data.rdm_data.y_axis = dummy_y_rdm_axis;
    % Plotly-specific properties for heatmap (optional, can be handled in JS layout)
    % These are hints for Plotly if you directly pass this as a trace in MATLAB
    % but for JS, you mostly need the data structure.
    % rdm_data.type = 'heatmap';
    % rdm_data.colorscale = 'Viridis';
    % rdm_data.colorbar = struct('title', 'Intensity');
    output_graph_data{end+1} = rdm_data;

    % 2. Dummy JRC Raw Breath (Scatter plot)
    % Expecting 'x_data' and 'y_data'
    dummy_raw_breath_y = sin(dummy_time * 2 * pi * 0.2) + randn(num_points, 1) * 0.1; % Ensure column vector
    raw_breath_data = struct();
    raw_breath_data.name = 'JRC Raw Breath';
    raw_breath_data.x_data = dummy_time;
    raw_breath_data.y_data = dummy_raw_breath_y;
    raw_breath_data.type = 'scatter';
    raw_breath_data.mode = 'lines';
    raw_breath_data.line = struct('color', '#1f77b4');
    output_graph_data{end+1} = raw_breath_data;

    % 3. Dummy JRC Filtered Breath (Scatter plot)
    % Expecting 'x_data' and 'y_data'
    dummy_filtered_breath_y = sin(dummy_time * 2 * pi * 0.2) * 0.8; % Smoother sine wave
    filtered_breath_data = struct();
    filtered_breath_data.name = 'JRC Filtered Breath';
    filtered_breath_data.x_data = dummy_time;
    filtered_breath_data.y_data = dummy_filtered_breath_y;
    filtered_breath_data.type = 'scatter';
    filtered_breath_data.mode = 'lines';
    filtered_breath_data.line = struct('color', '#ff7f0e');
    output_graph_data{end+1} = filtered_breath_data;

    % 4. Dummy JRC Breath PSD (Scatter plot - Frequency vs. Power)
    % Expecting 'x_data' and 'y_data'
    dummy_psd_breath_y = exp(-(dummy_freq_psd - 0.2).^2 / 0.05) + 0.1 * randn(num_points, 1); % Add noise and peak
    breath_psd_data = struct();
    breath_psd_data.name = 'JRC Breath PSD';
    breath_psd_data.x_data = dummy_freq_psd;
    breath_psd_data.y_data = dummy_psd_breath_y;
    breath_psd_data.type = 'scatter';
    breath_psd_data.mode = 'lines';
    breath_psd_data.line = struct('color', '#2ca02c');
    output_graph_data{end+1} = breath_psd_data;

    % 5. Dummy JRC Heartbeat PSD (Conditional based on hearbeat_test)
    % Expecting 'x_data' and 'y_data'
    if hearbeat_test
        dummy_psd_heart_y = exp(-(dummy_freq_psd - 1.2).^2 / 0.05) + 0.05 * randn(num_points, 1); % Add noise and peak
        heart_psd_data = struct();
        heart_psd_data.name = 'JRC Heartbeat PSD';
        heart_psd_data.x_data = dummy_freq_psd;
        heart_psd_data.y_data = dummy_psd_heart_y;
        heart_psd_data.type = 'scatter';
        heart_psd_data.mode = 'lines';
        heart_psd_data.line = struct('color', '#9467bd');
        output_graph_data{end+1} = heart_psd_data;
    end

    % 6. Dummy Sensor Breath PSD (Conditional based on sensor_test)
    % Expecting 'x_data' and 'y_data'
    if sensor_test
        dummy_psd_sensor_y = exp(-(dummy_freq_psd - 0.3).^2 / 0.05) + 0.1 * randn(num_points, 1);
        sensor_psd_data = struct();
        sensor_psd_data.name = 'Sensor Breath PSD';
        sensor_psd_data.x_data = dummy_freq_psd;
        sensor_psd_data.y_data = dummy_psd_sensor_y;
        sensor_psd_data.type = 'scatter';
        sensor_psd_data.mode = 'lines';
        sensor_psd_data.line = struct('color', '#8c564b');
        output_graph_data{end+1} = sensor_psd_data;

        % 7. Dummy Respirator Monitor (Time vs. Temp)
        % Expecting 'x_data' and 'y_data'
        dummy_respirator_temp_y = 25 + sin(dummy_time * 2 * pi * 0.1) * 2 + randn(num_points, 1) * 0.2;
        respirator_monitor_data = struct();
        respirator_monitor_data.name = 'Respirator Monitor';
        respirator_monitor_data.x_data = dummy_time;
        respirator_monitor_data.y_data = dummy_respirator_temp_y;
        respirator_monitor_data.type = 'scatter';
        respirator_monitor_data.mode = 'lines';
        respirator_monitor_data.line = struct('color', '#e377c2');
        output_graph_data{end+1} = respirator_monitor_data;
    end


    % --- 3. Generate Dummy Analysis Results (scalar values) ---
    output_analysis_results.jrcBreathRate = 12.5 + rand() * 2; % breaths/min
    output_analysis_results.range = 3.2 + rand() * 0.5; % meters
    output_analysis_results.speed = 0.15 + rand() * 0.05; % m/s

    if hearbeat_test
        output_analysis_results.jrcHeartRate = 70 + rand() * 5; % beats/min
    end
    if sensor_test
        output_analysis_results.sensorBreathRate = 11.8 + rand() * 1.5; % breaths/min
    end

    fprintf('Dummy data generation complete.\n');

catch ME
    % --- 4. Error Handling ---
    % Capture any errors during execution.
    errors = ME.message;
    % Get a more detailed report for better debugging:
    errors = [errors, newline, getReport(ME, 'extended', 'hyperlinks', 'off')];
    fprintf('Error in runAnalysisWrapper (DUMMY MODE): %s\n', errors);
    % Ensure outputs are empty on error
    output_graph_data = {};
    output_analysis_results = struct();
end

end