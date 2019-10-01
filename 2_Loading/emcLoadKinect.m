function [ tsvFile ] = emcLoadKinect( filename, cfg )
% Loads a single tsv file from the filename into a structure
% 
% syntax
% tsvFile = emcLoadKinect(filename, cfg);
% 
% input parameters
% filename: string containing the path or filename of the file to load
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.makersToKeep: cell array containing the markers names to be kept in the tsv structure
%     *.fillgapFlag: boolean indicating if use fillgap or not (from Motion Capture Toolbox)
%     *.classificationFlag: boolean indicating if emcClassification should be run
%     *.absolutePositionFlag: boolean indicating if emcAbsolutePosition should be run
%     *.removePrefixFlag: boolean indicating if emcRemovePrefix should be run
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

%% COMPUTATION AREA
data = readtext(filename, ' ');
data(1:3,:) = []; %% Useless lines
numCol = size(data,2);
numLine = size(data,1);
numMarker = (numCol-2)/8;

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
tsvFile.markerName = cell(numMarker,1);
tsvFile.data = zeros(numLine,numMarker*3);
tsvFile.analogdata = [];
tsvFile.other.descr = 'DESCRIPTION	--';
tsvFile.other.timeStamp = 'TIME_STAMP	--';
tsvFile.other.dataIncluded = '3D';
tsvFile.info.filename = filename;

% Assign data
time = zeros(numLine,1);
for lineIdx = 1:numLine
    time(lineIdx) = data{lineIdx,1};
    markerIdx = 1;
    for colIdx = 3:8:numCol
        if lineIdx == 1
            tsvFile.markerName{markerIdx} = data{lineIdx,colIdx};
        end
        startMarker = colIdx + 2; % Skip the markername and 'point'
        endMarker = startMarker + 2; % 3 dimensions
        idx = 1;
        for i = startMarker:endMarker
            tsvFile.data(lineIdx, (markerIdx-1)*3+idx) = data{lineIdx,i};
            idx = idx + 1;
        end
        markerIdx = markerIdx + 1;
    end
end

% Calculate Frequency
diff = time(2:end)-time(1:end-1);
tsvFile.freq = round(1/mean(diff));
end

