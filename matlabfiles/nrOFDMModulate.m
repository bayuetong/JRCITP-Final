% nrOFDMModulate.m
function waveform = nrOFDMModulate(carrier, grid)
    waveform = complex(randn(carrier.SampleRate / 100, 1), randn(carrier.SampleRate / 100, 1)); % Dummy short waveform
end