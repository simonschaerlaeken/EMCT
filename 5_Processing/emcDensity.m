function [ tsvFile ] = emcDensity(tsvFile, cfg)
% Computes the density - the average distance of all markers from the
% baricenter
% 
% syntax
% tsvFile = emcDensity(tsvFile, cfg)
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
% tsvFile = emcDensity(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.density
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
disp('[PROCESSING] Density')
% Select markers
featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featMarker);
tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
% Find Centroid
centroidData = centroid3D(tsvFileMarker.data);
% Calculate distance from that centroid
density = zeros(size(tsvFileMarker.data,1),tsvFileMarker.nMarkers);
for j = 1:tsvFileMarker.nMarkers
    %tsv_file_marker = mcgetmarker(tsv_file,j);
    indexStart = ((j-1)*3)+1;
    indexEnd = indexStart+2;
    for i = 1:size(tsvFileMarker.data,1)
        density(i,j) = distancePoints3d(tsvFileMarker.data(i,indexStart:indexEnd),centroidData(i,:));
    end
end
% Average
density = nanmean(density,2);
if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Density'],'NumberTitle','off')
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        plot(time, density)
    else
        plot(density)
    end
end

% Return the final value
tsvFile.processing.density = density;
end


