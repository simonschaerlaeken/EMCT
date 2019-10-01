function [ tsvFile ] = emcCumulativeDistance(tsvFile, cfg)
% Computes the cumulative distance
% 
% syntax
% tsvFile = emcCumulativeDistance(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.featMarker: cell array containing the markernames on which the
%     feature is calculated (default: all)
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.featMarker = {'a','b','c'};
% cfg.display = true
% cfg.displayUnit = 's'
% tsvFile = emcCumulativeDistance(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.cumuldist
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
% featMarker - markers on which the feature is calculated
if ~isfield(cfg, 'featMarker')||isempty(cfg.featMarker) 
    disp('cfg.featMarker was not setup. Default: all markers')
    cfg.featMarker = tsvFile.markerName;
end
% display - if true, plot the result
if ~isfield(cfg, 'display')
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
end
% displayUnit - str containing the unit of duration - 's' or 'tf'
if ~isfield(cfg, 'displayUnit')
    disp('cfg.displayUnit was not setup. Default: "tf"')
    cfg.displayUnit = 'tf';
end

%% COMPUTATION AREA
disp('[PROCESSING] Cumulative Distance')
% Select markers
featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featMarker);
tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
tsvFileDist = mccumdist(tsvFileMarker);
% Display
if cfg.display
    colorset = distinguishable_colors(length(cfg.featMarker));
    figure('Name',[tsvFile.info.filename, ' - Cumulative Distance'],'NumberTitle','off')
    hold on
    for i = 1:size(tsvFileDist.data,2)
        if strcmp(cfg.displayUnit, 's')
            time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
            plot(time, tsvFileDist.data(:,i), 'color', colorset(i,:))
        else
            plot(tsvFileDist.data(:,i), 'color', colorset(i,:))
        end
    end
    hold off
    legend(cfg.featMarker)
end
% Return the final value
tsvFile.processing.cumuldist = tsvFileDist.data;
end


