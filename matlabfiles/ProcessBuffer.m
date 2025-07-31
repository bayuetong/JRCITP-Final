% ProcessBuffer.m
function [sensor_reading, Pt_start_Buf, Pt_stop_Buf] = ProcessBuffer(sensor_raw, Pt_start_Buf, Pt_stop_Buf)
    sensor_reading = rand() * 100; % Dummy sensor reading
    Pt_start_Buf = 1;
    Pt_stop_Buf = 15;
end