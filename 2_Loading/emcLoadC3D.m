function [tsvFile] = emcLoadC3D( filename, cfg )
% Opens a C3D file and create a tsvFile Structure
% afterwards
% 
% syntax
% [tsvFile] = emcLoadC3D( filename, cfg )
% 
% input parameters
% filename: string containing the path or filename of the file to convert
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL] 
%     *.deleteMarker: cell array of markers names to be deleted
%     *.changeMarker: cell array of markers names to be changed
%     *.fillgapFlag: boolean indicating if use fillgap or not (from Motion Capture Toolbox)
%   
% output
% []
% 
% examples
% cfg.fillgapFlag = true;
% [tsvFile] = emcLoadC3D( filename, cfg )
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
% Location check
if ~exist(filename, 'file')
    error('File does not exist')
elseif isdir(filename) % Directory
    error('Folder is not yet taken into account in this version of the toolbox/ Please use file path as input.')
end

% deleteMarker: cell array of markers names to be deleted
if ~isfield(cfg, 'deleteMarker')
    cfg.deleteMarker = {};
elseif ischar(cfg.deleteMarker)
    cfg.removeMarker = {cfg.deleteMarker};
end

% fillgapFlag: boolean indicating if use fillgap or not
if ~isfield(cfg, 'fillGapFlag')
    cfg.fillGapFlag = false;
end

%% COMPUTING AREA

tsvFile = mcread(filename);


% Delete marker
if ~isempty(cfg.deleteMarker)
    % Find the position of the mispelled string
    indexDeleteList = findIndexList(tsvFile.markerName, cfg.deleteMarker);
    
    % Create vector with all markers excepts those ones
    indexList = 1:numel(tsvFile.markerName);
    tmp = ones(1,numel(tsvFile.markerName));
    tmp(indexDeleteList) = 0;
    indexList = indexList(find(tmp));
    
    % Get markers
    tsvFile = mcgetmarker(tsvFile, indexList);
end

% Fill gaps
if cfg.fillGapFlag
    tsvFile = mcfillgaps(tsvFile, 'fillall');
end


end

