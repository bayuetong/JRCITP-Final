% dummy_img_propagation.m
function [carrier, WaveImg1, dmrsAntSymbols, dmrsAntIndices, dmrsSymbols_all, dmrsIndices_all, trBlk_tx, trBlkSizes, pdsch_all, pdschIndices, pdschIndices_all, encodeDLSCH] = dummy_img_propagation(img_padded_length, img_input, NSlots, carrier, simLocal, pdsch, encodeDLSCH, decodeDLSCHLocal, pdschextra, harqEntity)
    disp('Running dummy img_propagation...');
    waveform = [];
    dmrsSymbols_all = [];
    dmrsIndices_all = [];
    trBlk_tx = [];
    pdsch_all = cell(1, NSlots);
    pdschIndices_all = [];

    for nslot = 0:NSlots-1
        carrier.NSlot = nslot;
        [pdschIndices, pdschIndicesInfo] = nrPDSCHIndices(carrier, pdsch);
        pdschIndices_all = [pdschIndices_all pdschIndices];
        trBlkSizes = nrTBS(pdsch.Modulation, pdsch.NumLayers, numel(pdsch.PRBSet), pdschIndicesInfo.NREPerPRB, pdschextra.TargetCodeRate, pdschextra.XOverhead);

        trBlk = img_input(:, nslot + 1);
        trBlk_tx = [trBlk_tx trBlk];
        setTransportBlock(encodeDLSCH, trBlk, 0, harqEntity.HARQProcessID); % Assuming single codeword

        codedTrBlocks = encodeDLSCH(pdsch.Modulation, pdsch.NumLayers, pdschIndicesInfo.G, harqEntity.RedundancyVersion, harqEntity.HARQProcessID);

        pdschGrid = nrResourceGrid(carrier, simLocal.NTxAnts);
        pdschSymbols = nrPDSCH(carrier, pdsch, codedTrBlocks);
        pdschGrid(pdschIndices) = pdschSymbols;

        dmrsSymbols = nrPDSCHDMRS(carrier, pdsch);
        dmrsIndices = nrPDSCHDMRSIndices(carrier, pdsch);
        dmrsAntSymbols = dmrsSymbols;
        dmrsAntIndices = dmrsIndices;
        pdschGrid(dmrsAntIndices) = dmrsAntSymbols;
        dmrsSymbols_all = [dmrsSymbols_all dmrsSymbols];
        dmrsIndices_all = [dmrsIndices_all dmrsIndices];

        txWaveform = nrOFDMModulate(carrier, pdschGrid);
        waveform = [waveform; txWaveform];
        pdsch_all{nslot+1} = pdsch; % Store pdsch config for this slot
    end

    WaveImg1 = resample(waveform, 2000000, 1228800); % Dummy resample
    % save '~/Desktop/rfnoc_test/U1Wave_Image_QPSK50M_200Msps.mat' WaveImg1; % Don't actually save

    disp('Dummy propagating image bitstream...');
    % system('python3 ~/Desktop/rfnoc_test/image_rfnoc_replay.py'); % Don't call python
end