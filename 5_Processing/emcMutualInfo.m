function [ tsvFile ] = emcMutualInfo( tsvFile, cfg)
%MUTUALINFO Function that return the mutual infro between the columns of
%a data matrix
%
% syntax
% mutualInfo = emcMutualInfo(data, cfg)
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
% mutualInfo = emcMutualInfo(data, cfg);
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
if ~isfield(cfg, 'title')
    disp('cfg.title was not setup.')
    cfg.title = 'Mutual Info';
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
disp('[PROCESSING] Mutual Info')
data = tsvFile.data;

nbCol = size(data,2);
mutualInfo = zeros(nbCol);

% Define the permutations between all the different columns
permutationMarker = nchoosek(1:nbCol,2);
% Computation of function for all the pairs
for pairIdx = 1:size(permutationMarker,1)
    % Select two columns
    idx1 = permutationMarker(pairIdx,1);
    idx2 = permutationMarker(pairIdx,2);
    mi = MutualInformation(data(:,idx1), data(:,idx2));
    % Average
    mutualInfo(idx1,idx2) = mi;
    mutualInfo(idx2,idx1) = mi;
end
% Computation of function for the identity
for i = 1:nbCol
    mi = MutualInformation(data(:,i), data(:,i));
    mutualInfo(i,i) = mi;
end

if cfg.display
    figure('Name',['Mutual information - ', cfg.title],'NumberTitle','off')
    imagesc(mutualInfo)
    set(gca,'XTick',1:length(cfg.colNames))
    set(gca,'YTick',1:length(cfg.colNames))
    set(gca,'xticklabel',cfg.colNames)
    set(gca,'yticklabel',cfg.colNames)
    rotateticklabel_imagesc(gca,270);
    colorbar
end
% Return the final value
tsvFile.processing.mutualInfo = mutualInfo;
end

