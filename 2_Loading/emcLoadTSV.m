function [ tsvFile ] = emcLoadTSV(filename, cfg)
% Loads a single tsv file from the filename into a structure
% 
% syntax
% tsvFile = emcLoadTSV(filename, cfg);
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
% tsvFile1 = emcLoadTSV('motion1.tsv', cfg);
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
% tsvFile2 = emcLoadTSV('motion2.tsv', cfg);
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
% fillGapFlag - If TRUE, will fill the gap in the motion capture file by
% extrapoling values before and after the gap
if ~isfield(cfg,'fillGapFlag')
    disp('INFO: fillGapFlag set up to TRUE')
    cfg.fillGapFlag = true;
end
% classificationFlag - if TRUE, classification will be made
if ~isfield(cfg,'classificationFlag')
    disp('INFO: No classification will be done')
    cfg.classificationFlag = false;
elseif ~islogical(cfg.classificationFlag)
    error('cfg.classificationFlag must be boolean')
end
% absolutePositionFlag - if TRUE, store the absolute position of each body
if ~isfield(cfg,'absolutePositionFlag')
    disp('INFO: Absolute value not stored')
    cfg.absolutePositionFlag = false;
elseif ~islogical(cfg.absolutePositionFlag)
    error('cfg.absolutePositionFlag must be boolean')
end
% makersToKeep - List of Markers that will be kept during the import. If
% list is empty, all markers will be kept
if ~isfield(cfg,'makersToKeep')
    disp('INFO: All markers will be kept')
    cfg.markersToKeep = {};
end
% removePrefix - str contraining the prefix to be removed 
if ~isfield(cfg,'removePrefixFlag')
    cfg.removePrefixFlag = false;
elseif ~islogical(cfg.removePrefixFlag)
    error('cfg.removePrefixFlag must be boolean')
end

% Location check
if ~exist(filename, 'file')
    error('File does not exist')
elseif isdir(filename) % Directory
    error('Folder is not yet taken into account in this version of the toolbox/ Please use file path as input.')
end


%% COMPUTATION AREA
%--------------------------------------------------------------------------
% load data (i.e the tsv files)
%--------------------------------------------------------------------------
% Load data
tsvFile = mcread(filename);
% Remove Prefix
if cfg.removePrefixFlag
    tsvFile = emcRemovePrefix(tsvFile, cfg);
end

% Keep markers we want
if ~isempty(cfg.markersToKeep)
    markersToKeepIdx = findindexlist(tsvFile.markerName, cfg.makersToKeep);
    tsvFile = mcgetmarker(tsvFile, markersToKeepIdx);
end

% Fill Gaps
if cfg.fillGapFlag
    tsvFile = mcfillgaps(tsvFile,'fillall');
end

%% INFOS
% Retrieve Classification infos
if cfg.classificationFlag
    tsvFile = emcClassification(tsvFile, cfg);
end

% Absolute Position
if cfg.absolutePositionFlag
    tsvFile = emcAbsolutePosition(tsvFile, cfg);
end

% Filename
filename = strsplit(tsvFile.filename, filesep);
filename = filename{end};
tsvFile.info.filename = filename;

end

