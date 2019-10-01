function [] = emcSaveContinuousFeature( tsvFile, cfg )
%EMCSAVECONTINUOUSFEATURE Summary of this function goes here
%   Detailed explanation goes here

time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
matrix = time'; 
headers = ['Time', cfg.exportContinuousFeatureList]';

for feat_idx = 1:numel(cfg.exportContinuousFeatureList)
    feature_name = cfg.exportContinuousFeatureList{feat_idx};
    matrix = horzcat(matrix,mean(tsvFile.processing.(feature_name),2));
end

csvwrite_with_headers([cfg.outputDir, filesep,'FeaturesMocap_', tsvFile.info.filename,'.csv'],matrix,headers)
end

