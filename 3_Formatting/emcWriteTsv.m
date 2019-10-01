function [] = emcWriteTsv( tsvFile, cfg )
% Write a TSV file and can format
% 
% syntax
% emcDivideMultipleBodies(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     * -
%     [OPTIONAL]
%     *.outputDir: string containing the path to the output directory
%     *.orderMarkerFlag: boolean saying if emcOrderMarker should be
%     computed
%     *.removePrefixFlag: boolean saying if emcRemovePrefix should be
%     computed
%     *.changeMarkerFlag: boolean saying if emcChangeMarkerName should be
%     computed 
%     *.deleteMarker: cell array of markers to be deleted
%     *.fillgapFlag: boolean indicating if use fillgap or not
% output
% []
% 
% examples
% cfg.outputDir = 'C:\Users\Desktop';
% cfg.removePrefixFlag = true;
% cfg.removePrefix = {'BODY1'}
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

% removePrefixFlag: boolean saying if emcRemovePrefix should be computed
if ~isfield(cfg,'removePrefixFlag')
    cfg.removePrefixFlag = false;
elseif ~islogical(cfg.removePrefixFlag)
    error('cfg.removePrefixFlag must be boolean')
elseif cfg.removePrefixFlag % if true
    if ~isfield(cfg,'removePrefix')
        error('cfg.removePrefix must be define')
    end
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

% fillgapFlag: boolean indicating if use fillgap or not
if ~isfield(cfg, 'fillGapFlag')
    cfg.fillGapFlag = false;
end

%% COMPUTATION AREA
% Remove Prefix
if cfg.removePrefixFlag
    cfg.removePrefix = {cfg.removePrefix};
    tsvFile = emcRemovePrefix(tsvFile, cfg);
end

% Change Marker
if cfg.changeMarkerFlag
    tsvFile = emcChangeMarkerName(tsvFile, cfg);
end

% Change order marker
if cfg.orderMarkerFlag
    tsvFile = emcOrderMarker(tsvFile, cfg);
end
% Change order marker
if cfg.fillGapFlag
    tsvFile = mcfillgaps(tsvFile, 'fillall');
end
% Write TSV
mcwritetsv(tsvFile,cfg.outputDir)


