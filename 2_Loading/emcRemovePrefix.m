function [ tsvFile ] = emcRemovePrefix( tsvFile, cfg )
% Removes a prefix / different prefixes to the list of markers of a tsvFile
% 
% syntax
% tsvFile = emcRemovePrefix(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.removePrefix: string or cell array containing the prefix
%     [OPTIONAL]
%     -
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.removePrefix = {'body1_'];
% tsvFile = emcRemovePrefix(tsvFile, cfg);
% 
% comments
% stores the prefix in tsvFile.info.removePrefix
% 
% see also
% emcLoadSingle
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
errorIfNotField(cfg, 'removePrefix')
% Prepare the cfg.removePrefix as a list
if iscell(cfg.removePrefix)
    % GOOD
elseif isstring(cfg.removePrefix)
    cfg.removePrefix = {cfg.removePrefix};
else
    error('Prefix should be either a string or a cell array containing strings')
end
%% COMPUTING AREA
% Compute for every cfg.removePrefix
disp('[INFO] Removing prefix')
for prefixIdx = 1:numel(cfg.removePrefix)
    prefixTmp = cfg.removePrefix{prefixIdx};
    [tsvFile.markerName] = removePrefixList(tsvFile.markerName, prefixTmp);
    tsvFile.info.removePrefix = prefixTmp;
end

