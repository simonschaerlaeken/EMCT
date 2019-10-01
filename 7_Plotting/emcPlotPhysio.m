function [] = emcPlotPhysio( tsvFile, cfg )
% Plots the physio data contained in a tsvFile
% 
% syntax
% emcPlotPhysio(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.physioList: cell array containing the feature names to be plotted
%     (default all)
%     *.featureName: string containing the name of the feature you might
%     want to plot
%     
% output
% -
% 
% examples
% cfg.physioList = {'ACC','BVP','EDA','HR'};
% emcPlotPhysio(tsvFile, cfg);
% 
% comments
% -
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
% filename - for title
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end

% Check if physio exists
if ~isfield(tsvFile, 'analogdata') || isempty(tsvFile.analogdata)
    error('Physiological data must be added before plotting the physio')
end

% physioList - cell array containing the physio names to be plotted
if ~isfield(cfg, 'physioList')
   cfg.physioList = {tsvFile.analogdata(:).physioName};
elseif ~iscell(cfg.physioList)
   cfg.physioList = {cfg.physioList};
end
% Check physioList in the phyio
for physioIdx = 1:numel(cfg.physioList)
    physioName = cfg.physioList{physioIdx};
    if ~ismember({tsvFile.analogdata(:).physioName},physioName)
        error(['Physio data ', physioName, 'does not exist in tsvFile.analogdata'])
    end
end

% Check if feature exists
if isfield(cfg, 'featureName') && ~isempty(cfg.featureName)
    if ~isfield(tsvFile, 'processing')
        error('Processing must be computed before plotting features')
    end
    if ~isfield(tsvFile.processing, cfg.featureName)
        error(['Feature ', cfg.featureName, 'does not exist in tsvFile.processing'])
    end
else
    cfg.featureName = '';
end

%% WORKING AREA
colorset = distinguishable_colors(length(cfg.physioList));
%title = [tsvFile.info.filename, ' - ', derivName];
figure('Name','Physiological Data','NumberTitle','off')
hold on
legendName = {};
if ~isempty(cfg.featureName)
    ax1 = subplot(2,1,1);
    hold on
end
for physioIdx = 1:numel(cfg.physioList)
    % Feature Name
    physioName = cfg.physioList{physioIdx};
    % Find physio
    physioStructIdx = findIndexList({tsvFile.analogdata(:).physioName}, {physioName});
    % Feature Data
    physioData = tsvFile.analogdata(physioStructIdx).data; 
    % Plot
    for colIdx = 1:size(physioData,2)
        % frequency
        freq = tsvFile.analogdata(physioStructIdx).freq(colIdx);
        % time
        time = 0:(1/freq):((size(physioData,1)-1)/freq);
        % plot
        plot(time,physioData(:,colIdx), 'color', colorset(physioIdx,:));
        % legend
        if size(physioData,2) > 1
            legendName{end+1} = [physioName, num2str(colIdx)];
        else
            legendName{end+1} = physioName;
        end
    end
end
hold off
legend(legendName)

if ~isempty(cfg.featureName)
    ax2 = subplot(2,1,2);
    % Feature Data
    featureData = tsvFile.processing.(cfg.featureName);
    % Plot
    % one value
    if size(featureData,1) == 1
        disp(['Feature ', featureName, ' ='])
        disp(featureData)
    % vector
    else
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        plot(time, featureData);
    end
    if isfield(cfg, 'featMarker') && numel(cfg.featMarker) == size(featureData,2)
        legend(ax2, cfg.featMarker);
    elseif numel(tsvFile.markerName) == size(featureData,2)
        legend(ax2, tsvFile.markerName);
    end
end
end
