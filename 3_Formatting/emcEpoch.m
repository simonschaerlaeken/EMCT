function [ tsvFile ] = emcEpoch( tsvFile, cfg )
% Computes the derivative of the motion 
% 
% syntax
% tsvFileEpooh = emcEpoch(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.epochInOut: array of int indexes for the beginning of the sequence
%     and the end -- array must contain only 2 int: ex [2 150].
%     [OPTIONAL]
%     -
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.epochInOut = [2 150];
% tsvFile = emcEpoch(tsvFile, cfg);
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
errorIfNotField(cfg, 'epochInOut')
if numel(cfg.epochInOut) ~= 2
    error('Array must contain only 2 int')
elseif cfg.epochInOut(1) >= cfg.epochInOut(2)
    error('First value must be lower than second value')
elseif cfg.epochInOut(1) >= size(tsvFile.data,1) || cfg.epochInOut(2) > size(tsvFile.data,1)
    error('Out of boundaries')
end
%% COMPUTATION AREA
disp('[INFO] Epoch tsvFile')
tsvFile.data = tsvFile.data(cfg.epochInOut(1):cfg.epochInOut(2),:);
if isfield(tsvFile, 'rotationdata')
    tsvFile.rotationdata = tsvFile.rotationdata(cfg.epochInOut(1):cfg.epochInOut(2),:);
end
tsvFile.nFrames = size(tsvFile.data,1);
end

