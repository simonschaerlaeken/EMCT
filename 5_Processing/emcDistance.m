function [ tsvFile ] = emcDistance(tsvFile, cfg)
% Computes the distance between 2 markers or to the ground from one marker
% 
% syntax
% tsvFile = emcDistance(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.featDistMarker: cell array containing 2 markernames or 1
%     markernames and "ground"
%     [OPTIONAL]
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.featDistMarker = {'a','b'};
% cfg.display = true
% tsvFile = emcDistance(tsvFile, cfg);
% 
% cfg.featDistMarker = {'a','ground'};
% cfg.display = true
% cfg.displayUnit = 's'
% tsvFile = emcDistance(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.distance
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
errorIfNotField(cfg,'featDistMarker')
if numel(cfg.featDistMarker) ~= 2
    error('Please enter 2 markers in cfg.featMarker for distance to be computed between them')
else
    if strcmp(cfg.featDistMarker{1}, 'ground')
        % Calculate distance to ground
        groundFlag = true;
        cfg.featMarker = cfg.featDistMarker(2);
    elseif strcmp(cfg.featDistMarker{2}, 'ground')
        % Calculate distance to ground
        groundFlag = true;
        cfg.featMarker = cfg.featDistMarker(1);
    else
        % Calculate distance between markers
        groundFlag = false;
    end
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
if groundFlag
    disp('[PROCESSING] Distance from ground')
    % Find index marker
    featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featDistMarker);
    tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
    % Take distance from ground (zCoord)
    distance = tsvFileMarker.data(:,3);
    if cfg.display
        figure('Name',[tsvFile.info.filename, ' - Distance from ground'],'NumberTitle','off')
        plot(distance)
    end
else
    disp('[PROCESSING] Distance between markers')
    % Find index marker
    featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featDistMarker);
    tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
    % Calculate distance between markers
    distance = zeros(size(tsvFileMarker.data,1),1);
    %Compute distance for everu timeframe
    for i = 1:size(tsvFileMarker.data,1)
        distance(i) = distancePoints3d(tsvFileMarker.data(i,1:3),tsvFileMarker.data(i,4:6));
    end
    % Average
    distance = mean(distance,2);
    if cfg.display
        figure('Name',[tsvFile.info.filename, ' - Distance between markers'],'NumberTitle','off')
        if strcmp(cfg.displayUnit, 's')
            time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
            plot(time, distance)
        else
            plot(distance)
        end
    end
end

% Return the final value
tsvFile.processing.distance = distance;
end


