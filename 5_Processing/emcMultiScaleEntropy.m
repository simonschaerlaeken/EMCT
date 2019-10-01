function [ tsvFile ] = emcMultiScaleEntropy( tsvFile, cfg )
% Compute the multiscale entropy as defined by Costa 2003.
% 
% syntax
% tsvFile = emcMultiScaleEntropy(tsvFile,cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.scaleList: array of int of all the different scale for the multi scale entropy
%     *.featureName: feature to calculate the sample entrop on
%     [OPTIONAL]
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.scaleList = [1,2,3,4,5,6];
% cfg.featureName = 'speed';
% tsvFile = emcMultiScaleEntropy(tsvFile,cfg);
%
% comments
% the output is written in a field named mse
% the output is also written in a field named mseParam contianing the area
% below the curve and the coef of the first degree linear fit 
% 
% see also
% -
%
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
% scaleList: array of int of all the different scale for the multi scale entropy
errorIfNotField(cfg, 'scaleList')

% featureName: feature to calculate the sample entrop on
errorIfNotField(cfg, 'featureName');

% display - if true, plot the result
if ~isfield(cfg, 'display')
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
end
%% COMPUTATION AREA
% Create final table
mse = zeros(numel(cfg.scaleList),size(tsvFile.processing.(cfg.featureName),2));

for scaleIdx = 1:numel(cfg.scaleList)
    scale = cfg.scaleList(scaleIdx);
    for segmIdx=1:size(tsvFile.processing.(cfg.featureName), 2)
        vector=(tsvFile.processing.(cfg.featureName)(:,segmIdx));
        coarsegrainedVector = coarsegrain(vector, scale, 'mean');
        % define Sample Entropy Parameters
        dim=2;
        r=0.2*std(vector);
        tau=1;
        % compute SampEn for each column (= each marker)  
        mse(scaleIdx,segmIdx) = sampleEntropy(dim, r, coarsegrainedVector, tau);
    end
end
tsvFile.processing.mse = mse;
mseMean = nanmean(tsvFile.processing.mse,2);
if cfg.display
    figure('Name',[tsvFile.info.filename,' - ',cfg.featureName, ' - Multi Scale Sample Entropy'],'NumberTitle','off')
    plot(mseMean)
end
integr = trapz((1:numel(mseMean))',mseMean);
polynomeLin = polyfit((1:numel(mseMean))',mseMean,1);
tsvFile.processing.mseIntegral = integr;
tsvFile.processing.mseSlope = polynomeLin(1);
end

