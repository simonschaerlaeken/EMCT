function [ tsvFile ] = emcChangeMarkerName( tsvFile, cfg )
% Change the name of markers of a tsvFile
% 
% syntax
% tsvFile = emcChangeMarkerName(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.changeMarker: cell array of markers names to be changed
%     [OPTIONAL]
%     -
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.changeMarker = {'tutu','e'; 'baba','f'};
% tsvFile = emcChangeMarkerName(tsvFile, cfg);
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
% changeMarker: cell array of markers names to be changed
errorIfNotField(cfg, 'changeMarker')
if ~iscell(cfg.changeMarker) || size(cfg.changeMarker,2) ~= 2
    error('cfg.changeMarker must be formatted as a cell table of dimension 2 for columns and the number of markers to be changed in lines')
end
%% COMPUTING AREA
% Compute for every cfg.removePrefix
disp('[INFO] Changing markername')

% Change marker
for correctionIdx = 1:size(cfg.changeMarker,1)
    faultStr = cfg.changeMarker{correctionIdx,1};
    % Check if fault string is present
    if any(strcmp(tsvFile.markerName,faultStr))
        % Find the position of the mispelled string
        indexList = findIndexList(tsvFile.markerName, {faultStr});

        % Replace if found
        if numel(indexList)>1
            error('Multiple occurences of the same marker name');
        else
            index = indexList(1);
            if index ~= 0
                disp(['Replacing marker:' faultStr ' with ' cfg.changeMarker{correctionIdx,2}])
                tsvFile.markerName{index} = cfg.changeMarker{correctionIdx,2};
            end
        end
    end
end

end

