function [ tsvFile ] = emcLoadCSV( filename, cfg )
% Loads a single CSV file from the filename into a structure
% 
% syntax
% tsvFile = emcLoadCSV(filename, cfg);
% 
% input parameters
% filename: string containing the path or filename of the file to load
% cfg: configuration structure
%     [MANDATORY]
%     *.suffixPosition: string containing the suffix of the first position
%     markers
%     *.suffixRotation: string containing the suffix of the first rotation
%     markers
%     *.nbDimPosition: num representing the number of dimensions of position (usually 3)
%     *.nbDimRotation: num representing the number of dimensions of rotation (usually 4)
%     *.timeName: string containing the name of the header for the time
%     column
%     [OPTIONAL]
%     *.makersToKeep: cell array containing the markers names to be kept in the tsv structure
%     *.fillgapFlag: boolean indicating if use fillgap or not (from Motion Capture Toolbox)
%     *.classificationFlag: boolean indicating if emcClassification should be run
%     *.absolutePositionFlag: boolean indicating if emcAbsolutePosition should be run
%     *.removePrefixFlag: boolean indicating if emcRemovePrefix should be run
%     *.epoch: vector containing the start and end of epoch
%     *.rotationToAxisFlag: boolean indicating if the rotation should be
%     transform to axis
%     *.rotationType: str containing the rotation type ex "quaternion"
% output
% tsvFile: MoCap data structure
% 
% examples
% tsvFile1 = emcLoadKinect('motion1.tsv', cfg);
% 
% cfg.makersToKeep = {'a','b','c'};
% cfg.fillgapFlag = true;
% cfg.absolutePositionFlag = true;
% cfg.removePrefixFlag = true;
% cfg.removePrefix = {'body1_'};
% cfg.classificationFlag = true;
% cfg.classificationType = 'prefix';
% cfg.classInfo = {'className1','className2'};
% cfg.classPrefixPosition = {1,2}; 
% cfg.classListPrefix = {{'1','2','3'},...
%                        {'1','2'}};
% cfg.classListValue = {{'class1_1','class1_2','class1_3'},...
%                       {'class1_2','class2_2'}}; 
% tsvFile2 = emcLoadKinect('motion2.tsv', cfg);
% 
% comments
% -
% 
% see also
% emcClassification
% emcAbsolutePosition
% emcRemovePrefixFlag

% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
errorIfNotField(cfg, 'suffixPosition')
errorIfNotField(cfg, 'nbDimPosition')
errorIfNotField(cfg, 'timeName')

% Rotation
rotationFlag = true;
if ~isfield(cfg,'suffixRotation') || ~ischar(cfg.suffixRotation)
    warning('No rotation information will be retained')
    rotationFlag = false;
end
% epoch
cfg.epochFlag = false;
if isfield(cfg, 'epoch') && ~isempty(cfg.epoch)
    cfg.epochFlag = true;
end
% Rotation to axis
if isfield(cfg,'rotationToAxisFlag') && cfg.rotationToAxisFlag
    errorIfNotField(cfg, 'rotationType')
else
    cfg.rotationToAxisFlag = false;
end
%% COMPUTATION AREA
[dataNum, dataText] = swallow_csv(filename);
header = dataText(1,:); %% Useless lines
dataNum = dataNum(2:end,:);

% Assign Data
listIndexPosition = findIndexListEndsWith(header, cfg.suffixPosition);
numMarker = numel(listIndexPosition);
numCol = numMarker*cfg.nbDimPosition;
numLine = size(dataNum,1);
data = zeros(numLine,numCol);
markerName = cell(numMarker,1);
for i = 1:numMarker
    idxIn = ((i-1)*cfg.nbDimPosition)+1;
    idxOut = idxIn+cfg.nbDimPosition-1;
    data(:,idxIn:idxOut) = dataNum(:,listIndexPosition(i):listIndexPosition(i)+cfg.nbDimPosition-1);
    markerName{i} = header{listIndexPosition(i)}(1:end-length(cfg.suffixPosition));
end

% Assign Data Rotation
if rotationFlag
    listIndexPosition = findIndexListEndsWith(header, cfg.suffixRotation);
    numCol = numMarker*cfg.nbDimRotation;
    rotationdata = zeros(numLine,numCol);
    for i = 1:numMarker
        idxIn = ((i-1)*cfg.nbDimRotation)+1;
        idxOut = idxIn+cfg.nbDimRotation-1;
        rotationdata(:,idxIn:idxOut) = dataNum(:,listIndexPosition(i):listIndexPosition(i)+cfg.nbDimRotation-1);
    end
    rotationDim = cfg.nbDimRotation;
    if cfg.rotationToAxisFlag
        % Transform the rotation data into 3D axis
        if strcmp(cfg.rotationType, 'quaternion')
            % From quaternion
            xAxisOrig = [1 0 0];
            yAxisOrig = [0 1 0];
            zAxisOrig = [0 0 1];
            rotationdatatmp = zeros(size(rotationdata,1), numMarker*9);
            for markerIdx = 1:numMarker
                idxIn = ((markerIdx-1)*cfg.nbDimRotation)+1;
                idxOut = idxIn+cfg.nbDimRotation-1;
                idxRotIn = ((markerIdx-1)*9)+1; % 9 columns since 3 columns for each axis and 3 axes
                idxRotOut = idxRotIn+8;
                for i = 1:size(rotationdata,1)
                    xAxis = qRotatePoint(xAxisOrig,rotationdata(i,idxIn:idxOut)');
                    yAxis = qRotatePoint(yAxisOrig,rotationdata(i,idxIn:idxOut)');
                    zAxis = qRotatePoint(zAxisOrig,rotationdata(i,idxIn:idxOut)');
                    rotationdatatmp(i,idxRotIn:idxRotOut) = horzcat(xAxis',yAxis',zAxis');
                end
            end
            rotationdata = rotationdatatmp;
            rotationDim = 9;
        end
    end
end

% Create tsvFile
tsvFile.type = 'MoCap data';
tsvFile.filename = filename;
tsvFile.nFrames = numLine;
tsvFile.nCameras = 0;
tsvFile.nMarkers = numMarker;
tsvFile.freq = 0; % will be calculated later
tsvFile.nAnalog = 0;
tsvFile.anaFreq = 0;
tsvFile.timederOrder = 0;
tsvFile.markerName = markerName;
tsvFile.data = data;
if rotationFlag
    tsvFile.rotationdata = rotationdata; 
    tsvFile.rotationDim = rotationDim; 
    tsvFile.rotationMarkerName = markerName; 
end
tsvFile.analogdata = [];
tsvFile.other.descr = 'DESCRIPTION	--';
tsvFile.other.timeStamp = 'TIME_STAMP	--';
tsvFile.other.dataIncluded = '3D';
% Filename
filename = strsplit(tsvFile.filename, filesep);
filename = filename{end};
tsvFile.info.filename = filename;
% Calculate Frequency
indexList = findIndexList(header, {cfg.timeName});
time = dataNum(:,indexList(1));
diff = time(2:end)-time(1:end-1);
tsvFile.freq = round(1/mean(diff));
% Epoching
if cfg.epochFlag
    startEpoch = cfg.epoch(1);
    endEpoch = cfg.epoch(2);
%     if startEpoch<time(1)
%         disp('Smaller than the first time recorded')
%         startIdxEpoch = 1;
%     else
%         startIdxEpoch = find(time>startEpoch,1);
%     end
%     if endEpoch>time(end)
%         disp('Bigger than the last time recorded')
%         endIdxEpoch = length(time);
%     else
%         endIdxEpoch = find(time>endEpoch,1);
%     end
%     cfg.epochInOut = [startIdxEpoch endIdxEpoch];
    cfg.epochInOut = [startEpoch endEpoch];
    tsvFile = emcEpoch(tsvFile, cfg);
end
end

