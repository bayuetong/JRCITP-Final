% HARQEntity.m
classdef HARQEntity
    properties
        HARQProcessID
        RedundancyVersion
        NewData
        SequenceTimeout
        NumCodewords
    end
    methods
        function obj = HARQEntity(harqSequence, rvSeq, numCodewords)
            obj.HARQProcessID = harqSequence(1);
            obj.RedundancyVersion = rvSeq(1);
            obj.NewData = true(1, numCodewords);
            obj.SequenceTimeout = false(1, numCodewords);
            obj.NumCodewords = numCodewords;
        end
    end
end