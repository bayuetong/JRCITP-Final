% RunTest.m
function [M, N] = RunTest(strCOM)
    disp("Running dummy RunTest...");
    M = -1; % Dummy return for rfnoc_test
    N = rand(600, 1) * 0.1 + 0.5; % Dummy sensor data
    disp("Dummy Tests completed!");
end