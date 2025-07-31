% nrOFDMDemodulate.m
function grid = nrOFDMDemodulate(carrier, waveform)
    grid = complex(randn(carrier.NSizeGrid * 12, carrier.SymbolsPerSlot), randn(carrier.NSizeGrid * 12, carrier.SymbolsPerSlot));
end