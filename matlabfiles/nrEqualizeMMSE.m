% nrEqualizeMMSE.m
function [pdschEq, csi] = nrEqualizeMMSE(pdschRx, pdschHest, noiseEst)
    pdschEq = complex(randn(size(pdschRx)), randn(size(pdschRx)));
    csi = ones(size(pdschRx)); % Dummy CSI
end