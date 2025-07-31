% nrResourceGrid.m
function grid = nrResourceGrid(carrier, numAnts)
    grid = complex(randn(carrier.NSizeGrid * 12, carrier.SymbolsPerSlot, numAnts), randn(carrier.NSizeGrid * 12, carrier.SymbolsPerSlot, numAnts));
end