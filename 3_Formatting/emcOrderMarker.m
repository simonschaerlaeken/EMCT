function [ tsvFileOrdered ] = emcOrderMarker( tsvFile, cfg )
% Reorders tsvFile placing markers into a specific order
% 
% syntax
% tsvFileOrdered = emcOrderMarker( tsvFile, cfg );
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.orderMarker: cell array of markernames in the specific order
%     [OPTIONAL]
%     -
%     
% output
% tsvFileOrder: MoCap data structure (with marker in a different order)
% 
% examples
% cfg.orderMarker = {'b','a','c'};
% tsvFile = emcOrderMarker( tsvFile, cfg );
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
errorIfNotField(cfg, 'orderMarker');
if numel(tsvFile.markerName) ~= numel(cfg.orderMarker)
    warning('Size of markerName and ordered list does not match -- extra markers will be deleted')
end

%% COMPUTING AREA
disp('[INFO] Reorder markers')
indexListOrdered = findIndexList(tsvFile.markerName, cfg.orderMarker);
tsvFileOrdered = mcgetmarker(tsvFile, indexListOrdered);
end

