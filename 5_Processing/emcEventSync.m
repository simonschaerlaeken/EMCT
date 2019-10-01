function [ tsvFile] = emcEventSync( tsvFile, cfg )
%EVENTSYNC Function that return the synchronisation between the columns of
%a data matrix
%
% syntax
% tsvFile = emcEventSync(data, cfg)
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
% [tsvFile] = emcEventSync(tsvFile, cfg);
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
disp('[PROCESSING] Event Sync')
if strcmp(cfg.featureName, 'rawdata')
    data = tsvFile.data;
else
    data = tsvFile.processing.(cfg.featureName);
    
nbCol = size(data,2);
sync = zeros(nbCol);
delay = zeros(nbCol);
% Define the permutations between all the different columns
permutationMarker = nchoosek(1:nbCol,2);
% Computation of function for all the pairs
for pairIdx = 1:size(permutationMarker,1)
    % Select two columns
    idx1 = permutationMarker(pairIdx,1);
    idx2 = permutationMarker(pairIdx,2);
    [Q, q] = event_synchronization_on_peaks(data(:,idx1), data(:,idx2), 'tot', 1);
    % Average
    sync(idx1,idx2) = Q;
    sync(idx2,idx1) = Q;
    delay(idx1,idx2) = q;
    delay(idx2,idx1) = q;
end
% Computation of function for the identity
for i = 1:nbCol
    [Q, q] = event_synchronization_on_peaks(data(:,i), data(:,i), 'tot', 1);
    sync(i,i) = Q;
    delay(i,i) = q;
end

if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Sync'],'NumberTitle','off');
    imagesc(sync);
    set(gca,'XTick',1:length(cfg.colNames));
    set(gca,'YTick',1:length(cfg.colNames));
    set(gca,'xticklabel',cfg.colNames);
    set(gca,'yticklabel',cfg.colNames);
    rotateticklabel_imagesc(gca,270);
    colorbar;
    figure('Name',[tsvFile.info.filename, ' - Delay'],'NumberTitle','off');
    imagesc(delay);
    set(gca,'XTick',1:length(cfg.colNames));
    set(gca,'YTick',1:length(cfg.colNames));
    set(gca,'xticklabel',cfg.colNames);
    set(gca,'yticklabel',cfg.colNames);
    rotateticklabel_imagesc(gca,270);
    colorbar
end
% Return the final value
tsvFile.processing.sync = sync;
tsvFile.processing.delay = delay;
end
    
