function [ tsvFile ] = emcMutualInformation( tsvFile, cfg )
%EMCMUTUALINFORMATION Computes the mutual information of each column of a feature based on the equation
% of Claude Shannon
%
% syntax
% tsvFile = emcMutualInformation(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.featureName: str feature name on which the sample entropy will be
%     computed
%     [OPTIONAL]
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.featureName = 'kinEnerg';
% cfg.display = true;
% tsvFile = emcMutualInformation(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.entrop
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
% featureName - feature to calculate the sample entrop on
errorIfNotField(cfg, 'featureName');

% display - if true, plot the result
if ~isfield(cfg, 'display')
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
end


%% COMPUTATION AREA
disp('[PROCESSING] Entropy')
% Processing
nbCol = size(tsvFile.processing.(cfg.featureName),2);
mutualInfo = zeros(nbCol);
for col1 = 1:nbCol
    for col2 = 1:nbCol
        mutualInfo(col1,col2) = MutualInformation(tsvFile.processing.(cfg.featureName)(:,col1),tsvFile.processing.(cfg.featureName)(:,col2));
    end
end

% Display
if cfg.display
    figure('Name',['Mutual Information - ' cfg.featureName],'NumberTitle','off')               
    axisValues = fieldnames(tsvFile.markerName);
    imagesc(mutualInfo)
    set(gca,'XTick',1:length(axisValues))
    set(gca,'YTick',1:length(axisValues))
    set(gca,'xticklabel',axisValues)
    set(gca,'yticklabel',axisValues)
    rotateticklabel_imagesc(gca,270);
    colorbar
end
% Return the final value
tsvFile.processing.mutualInfo = mutualInfo;
end


