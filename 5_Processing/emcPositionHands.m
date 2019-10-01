function [ tsvFile ] = emcPositionHands(tsvFile, cfg)
% Computes the amount of time the hands spend in either the upper, middle or
% lower part of the body. The upper part corresponds to "above the head",
% the middle corresponds to "between the head and pelvis", lower part
% corresponds to "below the pelvis"
% 
% syntax
% tsvFile = emcPositionHands(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.headMarker: cell array containing one markername for the head
%     *.pelvisMarker: cell array containing one markername for the pelvis
%     *.handsMarker: cell array containing two markername for the two hands
%     [OPTIONAL]
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.headMarker = {'a'};
% cfg.pelvisMarker = {'a'};
% cfg.handsMarker = {'c','d'};
% cfg.display = true;
% cfg.displayUnit = 's'
% tsvFile = emcPositionHands(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.positionHands
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
% Check if not missing
errorIfNotField(cfg, 'headMarker');
errorIfNotField(cfg, 'pelvisMarker');
errorIfNotField(cfg, 'handsMarker');


%% COMPUTATION AREA
disp('[PROCESSING] Hands Position')

% Find index marker
featHeadMarkerIndex = findIndexList(tsvFile.markerName, cfg.headMarker);
featRootMarkerIndex = findIndexList(tsvFile.markerName, cfg.pelvisMarker);
featHandsMarkerIndex = findIndexList(tsvFile.markerName, cfg.handsMarker);
% % Extract markers from tsv
tsvFileHeadMarker = mcgetmarker(tsvFile,featHeadMarkerIndex);
tsvFileRootMarker = mcgetmarker(tsvFile,featRootMarkerIndex);
tsvFileHandsMarker = mcgetmarker(tsvFile,featHandsMarkerIndex);
% Extract data from tsv_markers
upperThreshold = tsvFileHeadMarker.data(:,3);
rootThreshold = tsvFileRootMarker.data(:,3);
middleThreshold = mean(horzcat(rootThreshold,upperThreshold),2);
handsElevation = tsvFileHandsMarker.data(:,[3 6]);
% check position hands:
upper = zeros(size(handsElevation,1),1);
middle = zeros(size(handsElevation,1),1);
lower = zeros(size(handsElevation,1),1);
for i = 1:size(handsElevation,1)
    for j = 1:size(handsElevation,2)
        if handsElevation(i,j) > upperThreshold
            upper(i) = upper(i) + 1;
        elseif handsElevation(i,j) > middleThreshold
            middle(i) = middle(i) + 1;
        else
            lower(i) = lower(i) + 1;
        end
    end
end
handsElevation = horzcat(upper,middle,lower);
tsvFile.processing.positionHands = handsElevation;
% Display
if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Hands Position'],'NumberTitle','off')
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        plot(time, handsElevation)
    else
        plot(handsElevation)
    end
    title('Hands Postion - 1:Upper, 2: Middle, 3: Lower')
end
end


