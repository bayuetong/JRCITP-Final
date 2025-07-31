function [analysisResults, plotData, imageOutputPath] = unified_usrp_analysis(runTests_in, sensor_test_in, hearbeat_test_in, speed_test_in, real_time_in, ts_in, img_filename_in, strCOM_in)
    % Assign input arguments to local variables
    runTests = runTests_in;
    sensor_test = sensor_test_in;
    hearbeat_test = hearbeat_test_in;
    speed_test = speed_test_in;
    real_time = real_time_in;
    ts = ts_in;
    img_filename = img_filename_in;
    strCOM = strCOM_in;

    % Initialize output structures
    analysisResults = struct();
    plotData = {};
    imageOutputPath = '';

    % --- Paste your entire unified_usrp_test_30k.m script content here ---
    % (Make sure all variables used in the script are either inputs, defined within the function, or loaded from accessible files.)

    % IMPORTANT: For animated plots (e.g., `h1`, `h2`, `h3` with `animatedline`),
    % instead of drawing live, collect all data points into arrays throughout the loop.
    % Example:
    % At the start of the function:
    % collected_jrc_breath_rate_time = [];
    % collected_jrc_breath_rate_val = [];
    % (and similar for sensor breath rate and JRC heart rate)
    % Inside your `for NUM=len+1:INUM` loop, replace `addpoints(h1,NUM*ts,p_rate_breath_JRC);` with:
    % collected_jrc_breath_rate_time(end+1) = NUM*ts;
    % collected_jrc_breath_rate_val(end+1) = p_rate_breath_JRC;
    % (Do this for all animated lines.)

    % --- End of your script content ---

    % 2. Capture Scalar Results:
    % Populate the 'analysisResults' struct with key numerical outcomes.
    analysisResults.jrc_range = range;
    analysisResults.jrc_speed = v_est;
    analysisResults.jrc_breath_rate = p_rate_breath_JRC;
    if hearbeat_test
        analysisResults.jrc_heart_rate = p_rate_heart_JRC;
    end
    if sensor_test
        analysisResults.sensor_breath_rate = p_rate_breath_Sensor;
        analysisResults.respirator_temp_raw = respiratorTemp; % Optional: raw sensor data
    end
    analysisResults.bit_error_rate = ber;

    % 3. Prepare Plot Data for Plotly.js:
    % Convert your MATLAB plot data into a structured format (e.g., cell array of structs)
    % that can be easily translated into Plotly.js graph objects by Python.
    % Each struct should contain `type`, `x`, `y`, `z` (for 3D plots), `title`, `xlabel`, `ylabel`.

    % Example for Range-Doppler Map:
    plotData{1} = struct('type', 'heatmap', 'x', v_array, 'y', r_array(1:20), 'z', abs(x_ar(1:20,:)).^2, ...
                         'title', 'JRC RDM', 'xlabel', 'Speed (m/s)', 'ylabel', 'Range (m)');

    % Example for JRC Raw Breath Signal:
    plotData{2} = struct('type', 'scatter', 'mode', 'lines', 'x', ts:ts:INUM*ts, 'y', x1, ...
                         'title', 'JRC Raw Breath Signal', 'xlabel', 'Time (s)', 'ylabel', 'Phase Change');

    % (Add structs for all other plots like JRC breath signal after filter, PSDs, and the collected time-series data for breath/heart rates.)

    % 4. Capture Image Output Path:
    % The script saves the reconstructed image. Return its path.
    % Ensure `img_extension` is available from your script's logic.
    imageOutputPath = ['~/Desktop/matlab/test_data/temp' img_extension];

end