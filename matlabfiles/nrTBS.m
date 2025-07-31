% nrTBS.m
function tbs = nrTBS(modulation, numLayers, numPRB, nREPerPRB, targetCodeRate, xOverhead)
    % A very simplified dummy TBS calculation
    bitsPerSymbol = 2; % QPSK
    if contains(modulation, '16QAM'), bitsPerSymbol = 4; end
    if contains(modulation, '64QAM'), bitsPerSymbol = 6; end
    if contains(modulation, '256QAM'), bitsPerSymbol = 8; end

    effectiveRE = numPRB * nREPerPRB - xOverhead;
    tbs = round(effectiveRE * bitsPerSymbol * targetCodeRate * numLayers);
    if tbs < 1, tbs = 100; end % Ensure a minimum size for dummy
end