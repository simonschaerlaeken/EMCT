function [ tsvFile ] = emcEntropy( tsvFile, cfg )
%EMCENTROPY Computes the entropy of a feature based on the equation
%of Claude Shannon
% 
% syntax
% tsvFile = emcEntropy(tsvFile, cfg)
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
% tsvFile = emcEntropy(tsvFile, cfg);
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
entrop = Entropy(tsvFile.processing.(cfg.featureName)); % Use toolbox Shannon from Will Dwinnell
% Display
if cfg.display
    figure('Name',[tsvFile.info.filename,' - ',cfg.featureName, ' - Sample Entropy'],'NumberTitle','off')
    % Set NaN to zeros
    entropZero = entrop;
    entropZero(isnan(entropZero)) = 0;
    plot(entropZero)
    if strcmp(cfg.featureName, 'kinEnerg')
        set(gca,'XTick',1:length(cfg.j2spar.segmentName))
        set(gca,'xticklabel',cfg.j2spar.segmentName)
    else
        set(gca,'XTick',1:length(tsvFile.markerName))
        set(gca,'xticklabel',tsvFile.markerName)
    end
    rotateticklabel2(gca,45);
end
% Return the final value
tsvFile.processing.entrop = entrop;
end
