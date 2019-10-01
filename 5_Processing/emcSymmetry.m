function [ tsvFile ] = emcSymmetry(tsvFile, cfg)
% Computes the symmetry between the left and right part of the body as
% defined by Piana 2013
% 
% syntax
% tsvFile = emcSymmetry(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.leftPrefix: string prefix of the left side of the body markernames
%     -- they should all start with this prefix
%     *.rightPrefix: string prefix of the right side of the body markernames
%     -- they should all start with this prefix
%     [OPTIONAL]
%     *.featMarker: cell array containing the markernames on which the
%     feature is calculated (default: all)
%     *.display: boolean deciding if a figure is to be plotted (default: true)
%     *.displayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.leftPrefix = 'L';
% cfg.rightPrefix = 'R';
% cfg.featMarker = {'Ra','Rb','Lc','Ld'};
% cfg.display = true
% cfg.displayUnit = 's'
% tsvFile = emcSymmetry(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.symmetry
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

errorIfNotField(cfg, 'leftPrefix');
errorIfNotField(cfg, 'rightPrefix');

%% COMPUTATION AREA
disp('[PROCESSING] Symmetry')
% Select markers
featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featMarker);
tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);

% Find Right and Left markers
leftMarkersList = {};
rightMarkersList = {};
for markerIdx = 1:length(tsvFileMarker.markerName)
    if length(tsvFileMarker.markerName{markerIdx}) < length(cfg.leftPrefix) || length(tsvFileMarker.markerName{markerIdx}) < length(cfg.rightPrefix)
        continue
    end
    % Left
    if strcmp(tsvFileMarker.markerName{markerIdx}(1:length(cfg.leftPrefix)),cfg.leftPrefix)
        leftMarkersList{end+1} = tsvFileMarker.markerName{markerIdx};
    % Right
    elseif strcmp(tsvFileMarker.markerName{markerIdx}(1:length(cfg.rightPrefix)),cfg.rightPrefix)
        rightMarkersList{end+1} = tsvFileMarker.markerName{markerIdx};
    end
end
% Verify markers
leftMarkersListFinal = {};
rightMarkersListFinal = {};
for markerLeftIdx = 1:length(leftMarkersList)
    markerLeft = leftMarkersList{markerLeftIdx};
    markerRight = markerLeft; % Should be similar but..
    markerRight = [cfg.rightPrefix,markerRight(length(cfg.leftPrefix)+1:end)]; % Put and R instead of L
    if sum(ismember(rightMarkersList, markerRight)) > 0
        leftMarkersListFinal{end+1} = markerLeft;
        rightMarkersListFinal{end+1} = markerRight;
    end
end
% Convert into indexes
leftMarkersIndex = findIndexList(tsvFileMarker.markerName,leftMarkersListFinal);
rightMarkersIndex = findIndexList(tsvFileMarker.markerName,rightMarkersListFinal);

% Compute Sym
symmetry = zeros(size(tsvFileMarker.data,1), length(rightMarkersListFinal));
for markerIdx = 1:length(rightMarkersIndex)
    rightTsvFile = mcgetmarker(tsvFile,rightMarkersIndex(markerIdx));
    leftTsvFile = mcgetmarker(tsvFile,leftMarkersIndex(markerIdx));
    centroidData = centroid3D(horzcat(rightTsvFile.data,leftTsvFile.data));
    numerator = (rightTsvFile.data-centroidData) - (leftTsvFile.data-centroidData);
    denominator = abs(rightTsvFile.data-centroidData) + abs(leftTsvFile.data-centroidData);
    for i = 1:size(rightTsvFile.data)
        symmetry(i,markerIdx) = numerator(i,:)/denominator(i,:);
    end    
end
% Average
symmetry = mean(symmetry,2);
if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Symmetry'],'NumberTitle','off')
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        plot(time, symmetry)
    else
        plot(symmetry)
    end
end
% Return the final value
tsvFile.processing.symmetry = symmetry;
end


