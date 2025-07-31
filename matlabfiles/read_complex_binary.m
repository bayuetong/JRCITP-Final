% read_complex_binary.m
function data = read_complex_binary(filename)
    % Return dummy complex data
    len = 1228800; % Default length used in the main script
    data = complex(randn(len, 1), randn(len, 1));
end