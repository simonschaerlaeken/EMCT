function [ tsvFile ] = emcSelect( tsvFile, cfg )
% Function to select part of the signal and analyse the physiological
% signal on those epoched parts
% 
% syntax
% [tsvFile] = emcSelect( tsvFile, cfg )
% 
% input parameters
% tsvFile: motion capture structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL] 
%     *.featureName: string containing the name of the feature you might
%     want to plot
%     *.epochLimit: (n x 2) table containing the in and out point of epoch
%     (if this argument is set, nEpoch should not be set)
%     *.nEpoch: int number of epochs (Mandatory if points not defined)
%
% output
% tsvFile: motion capture structure
% 
% examples
% cfg.featureName = 'speed'
% cfg.epochLimit = [34, 189; 289, 678];
% [tsvFile] = emcSelect( tsvFile, cfg )
%
% cfg.nEpoch = 3;
% [tsvFile] = emcSelect( tsvFile, cfg )
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
% filename - for title
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end

% Check if physio exists
if ~isfield(tsvFile, 'analogdata') || isempty(tsvFile.analogdata)
    error('Physiological data must be added before plotting the physio')
end
% Check if feature exists
if isfield(cfg, 'featureName') && ~isempty(cfg.featureName)
    if ~isfield(tsvFile, 'processing')
        error('Processing must be computed before plotting features')
    end
    if ~isfield(tsvFile.processing, cfg.featureName)
        error(['Feature ', cfg.featureName, 'does not exist in tsvFile.processing'])
    end
else
    cfg.featureName = '';
end

% Set limit or number of epoch
if ~isfield(cfg, 'nEpoch')
    errorIfNotField(cfg, 'epochLimit')
    if size(cfg.epochLimit,2)~=2
        error('Size of epochLimit must be n X 2')
    end
    cfg.selectFlag = false;
elseif ~isfield(cfg, 'epochLimit')
    errorIfNotField(cfg, 'nEpoch')
    cfg.selectFlag = true;
end

%% WORKING AREA
emcPlotPhysio(tsvFile, cfg) % Plot physio with eventual feature
if cfg.selectFlag
[X,~] = ginput(cfg.nEpoch*2);
cfg.epochLimit = (reshape(X,2,cfg.nEpoch))';
end

for physioIdx = 1:numel(tsvFile.analogdata)
    tsvFile.analogdata(physioIdx).epochdata = cell(size(cfg.epochLimit,1),1);
    tsvFile.analogdata(physioIdx).epochInfo = cell(size(cfg.epochLimit,1),1);
    tsvFile.analogdata(physioIdx).processing = struct();
    for epochIdx = 1:size(cfg.epochLimit,1)
        epochIn = round((cfg.epochLimit(epochIdx,1)*tsvFile.analogdata(physioIdx).freq(1))+1);
        if epochIn<1
            epochIn = 1;
        end
        epochOut = round((cfg.epochLimit(epochIdx,2)*tsvFile.analogdata(physioIdx).freq(1))+1);
        if epochOut > size(tsvFile.analogdata(physioIdx).data,1)
            epochOut = size(tsvFile.analogdata(physioIdx).data,1);
        end
        tsvFile.analogdata(physioIdx).epochInfo{epochIdx} = [epochIn,epochOut];
        tsvFile.analogdata(physioIdx).epochdata{epochIdx} = tsvFile.analogdata(physioIdx).data(epochIn:epochOut, :);
        tsvFile.analogdata(physioIdx).processing.mean = mean(tsvFile.analogdata(physioIdx).epochdata{epochIdx},1);
        tsvFile.analogdata(physioIdx).processing.std = std(tsvFile.analogdata(physioIdx).epochdata{epochIdx},1);
    end
end
end

