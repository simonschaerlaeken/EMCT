function [] = emcWriteTsvMultBodies( tsvFile, cfg )
% Divides a TSV file containing multiple bodies (represented by prefixe)
% into multiple files containing one body and save them
% 
% syntax
% emcDivideMultipleBodies(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.prefixBodyList: cell array of prefix name for every body
%     [OPTIONAL]
%     *.outputDir: string containing the path to the output directory
%     *.orderMarkerFlag: boolean saying if emcOrderMarker should be
%     computed
%     *.removePrefixFlag: boolean saying if emcRemovePrefix should be
%     computed
%     *.changeMarkerFlag: boolean saying if emcChangeMarkerName should be
%     computed 
%     *.deleteMarker: cell array of markers to be deleted
% output
% []
% 
% examples
% cfg.prefixBodyList = {'body1_','body2_'};
% cfg.outputDir = 'C:\Users\Desktop';
% cfg.removePrefixFlag = true;
% cfg.changeMarkerFlag = true;
% cfg.changeMarker = {'haha', 'b';'hehe','c'};
% cfg.orderMarkerFlag = true;
% cfg.orderMarker = ['LFHD','RFHD','LBHD','RBHB','CLAV','STRN','LSHO','RSHO','LELB','RELB','LWRA','RWRA','LASI','RASI','LKNE','RKNE','LANK','RANK'];
% emcDivideMultipleBodies(tsvFile, cfg );
% 
% comments
% Saves the different file as "filename_prefixBody.tsv"
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
if ~isfield(cfg, 'outputDir')
    disp('cfg.outputDir not define - default = current directory')
    cfg.outputDir = pwd;
elseif ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir)
end

% prefixBodyList - Contains every body prefix
errorIfNotField(cfg, 'prefixBodyList');

% removePrefixFlag: boolean saying if emcRemovePrefix should be computed
if ~isfield(cfg,'removePrefixFlag')
    cfg.removePrefixFlag = false;
elseif ~islogical(cfg.removePrefixFlag)
    error('cfg.removePrefixFlag must be boolean')
end

% orderMarkerFlag: boolean saying if emcOrderMarker should be computed
if ~isfield(cfg,'orderMarkerFlag')
    cfg.orderMarkerFlag = false;
elseif ~islogical(cfg.orderMarkerFlag)
    error('cfg.orderMarkerFlag must be boolean')
end

% changeMarkerFlag: boolean saying if emcChangeMarkerName should be computed
if ~isfield(cfg,'changeMarkerFlag')
    cfg.changeMarkerFlag = false;
elseif ~islogical(cfg.changeMarkerFlag)
    error('cfg.changeMarkerFlag must be boolean')
end


%% COMPUTATION AREA
% Divide for each body
for prefixIdx = 1:numel(cfg.prefixBodyList)
    prefixBody = cfg.prefixBodyList{prefixIdx};
    markerNameList = tsvFile.markerName;
    [indexList] = findIndexListStartsWith(markerNameList, prefixBody);
    tsvFileBody = mcgetmarker(tsvFile, indexList);
    
    % Change filename
    filename = strsplit(tsvFileBody.filename, '.tsv');
    prefixBodyName =regexprep(prefixBody, '[^\d\w~!@#$%^&()_\-{}.]*',''); % Make sure there is no forbidden character
    filename = [filename{1}, '_', prefixBodyName, '.tsv'];
    tsvFileBody.filename = filename;
    
    % Remove Prefix
    if cfg.removePrefixFlag
        cfg.removePrefix = {prefixBody};
        tsvFileBody = emcRemovePrefix(tsvFileBody, cfg);
    end
    
    % Change Marker
    if cfg.changeMarkerFlag
        tsvFileBody = emcChangeMarkerName(tsvFileBody, cfg);
    end
    
    % Change order marker
    if cfg.orderMarkerFlag
        tsvFileBody = emcOrderMarker(tsvFileBody, cfg);
    end
    
    % Write TSV
    mcwritetsv(tsvFileBody,cfg.outputDir)
end


