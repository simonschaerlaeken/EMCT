function [ tsvFile ] = emcCorrelation( tsvFile, cfg )
%EMCCORRELATION Function that return the correlation between the columns of
%a data matrix
%
% syntax
% [corr] = emcCorrelation(data, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.colNames: cell of str containing the name of the columns
%     *.title: str contraining the title of the analysis
%     *.display: boolean deciding if a figure is to be plotted (default: true)
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.colNames = tsvFile.markerNames;
% cfg.display = true;
% data = tsvFile.data;
% [corr] = emcCorrelation(data, cfg);
% 
% comments
% feature must be saved after function! Is not retain in the structure
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
% Initialisation
%% CHECKING AREA

% display - if true, plot the result
if ~isfield(cfg, 'display')
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
end

% title - if true, plot the result
if ~isfield(cfg, 'featureName')
    disp('cfg.featureName was not setup. Default: "rawdata"')
    cfg.featureName = 'rawdata';
end

% colNames - list of marker or feature used in analysis
if ~isfield(cfg, 'colNames')
    disp('cfg.colNames was not setup. Default: Numeric')
    cfg.colNames = cell(1,size(tsvFile.data,2));
    for i = 1:size(tsvFile.data,2)
        cfg.colNames{i} = num2str(i);
    end
end
%% COMPUTATION AREA
% Initialisation
disp('[PROCESSING] Correlation')
if strcmp(cfg.featureName, 'rawdata')
    data = tsvFile.data;
else
    data = tsvFile.processing.(cfg.featureName);

corr = corrcoef(data);

if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Correlation'],'NumberTitle','off')
    imagesc(corr)
    set(gca,'XTick',1:length(cfg.colNames))
    set(gca,'YTick',1:length(cfg.colNames))
    set(gca,'xticklabel',cfg.colNames)
    set(gca,'yticklabel',cfg.colNames)
    rotateticklabel_imagesc(gca,270);
    colorbar
end
% Return the final value
tsvFile.processing.corr = corr;
end

