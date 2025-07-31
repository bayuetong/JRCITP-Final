% nrPDSCH.m
function symbols = nrPDSCH(carrier, pdsch, codedTrBlocks)
    symbols = complex(randn(length(codedTrBlocks)/2, 1), randn(length(codedTrBlocks)/2, 1));
end