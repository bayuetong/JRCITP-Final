% nrExtractResources.m
function [pdschRx, pdschHest, varargout] = nrExtractResources(pdschIndices, rxGrid, estChannelGrid)
    pdschRx = rxGrid(pdschIndices);
    pdschHest = estChannelGrid(pdschIndices);
end