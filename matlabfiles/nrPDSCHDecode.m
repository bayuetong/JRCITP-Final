% nrPDSCHDecode.m
function [dlschLLRs, rxSymbols] = nrPDSCHDecode(carrier, pdsch, pdschEq, noiseEst)
    dlschLLRs = {complex(randn(length(pdschEq)*2, 1), randn(length(pdschEq)*2, 1))}; % Dummy LLRs (QPSK)
    rxSymbols = {pdschEq};
end