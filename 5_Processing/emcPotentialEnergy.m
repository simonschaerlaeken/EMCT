function [ tsvFile ] = emcPotentialEnergy(tsvFile, cfg)
% Computes the potential energy of the markers (based on mcpotenergy from
% the Motion Capture Toolbox, Copyright 2008, University of Jyvaskyla,
% Finland )
% 
% syntax
% tsvFile = emcPotentialEnergy(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.m2jpar: vector contraining the index of the markers for the
%     transformation to joints (see the Motion Capture Toolbox)
%     *.j2spar: vector contraining the index of the joints for the
%     transformation to segments (see the Motion Capture Toolbox)
%     *.mcspar: vector contraining the weights of all the segments (see the Motion Capture Toolbox) 
%     *.paramComputation: str flag defining the mean of concatenation of
%     the data - 'mean','sum','median'
%     [OPTIONAL]
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     *.featMarker: markers on which the feature is calculated
%
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.m2jpar.type = 'm2jpar';
% cfg.m2jpar.nMarkers = 22;
% cfg.m2jpar.markerNum = {[4,9,11,12],[4,11],[4,5,1],[1,2,3]};
% cfg.m2jpar.markerName = {'PELVIS', 'RPELV', 'RTHIGH', 'RTIB'};
% cfg.j2spar.type = 'j2spar';
% cfg.j2spar.rootMarker = 1;
% cfg.j2spar.frontalPlane = [1 2];
% cfg.j2spar.parent = [0 1 2 3];
% cfg.j2spar.segmentName = cfg.m2jpar.markerName;
% cfg.segmidx = [0 1 8 7];
% cfg.mcspar = mcgetsegmpar('Dempster', cfg.segmidx);
% cfg.paramComputation = 'mean';
% cfg.display = true;
% cfg.displayUnit = 's'
% tsvFile = emcPotentialEnergy(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.potEnerg
% 
% see also
% mcpotenergy
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
% display - if true, plot the result
if ~isfield(cfg, 'display')
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
end
% displayUnit - str containing the unit of duration - 's' or 'tf'
if ~isfield(cfg, 'displayUnit')
    disp('cfg.displayUnit was not setup. Default: "tf"')
    cfg.displayUnit = 'tf';
end
% featMarker - markers on which the feature is calculated
if ~isfield(cfg, 'featMarker')||isempty(cfg.featMarker) 
    disp('cfg.featMarker was not setup. Default: all markers')
    cfg.featMarker = tsvFile.markerName;
end

errorIfNotField(cfg, 'm2jpar') % Transform markers to joints
errorIfNotField(cfg, 'j2spar') % Transform joints to segments
errorIfNotField(cfg, 'mcspar') % Weights of segments
errorIfNotField(cfg, 'paramComputation') % how concatenate

%% COMPUTATION AREA
disp('[PROCESSING] Potential Energy')
% Create Joint and segment files
tsvFileJoint = mcm2j(tsvFile, cfg.m2jpar);
tsvFileSegment = mcj2s(tsvFileJoint, cfg.j2spar);
% Compute Potential energy
potEnrgVector = mcpotenergy(tsvFileJoint, tsvFileSegment, cfg.mcspar);

% Plot
if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Potential Energy'],'NumberTitle','off')
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        plot(time, potEnrgVector)
    else
        plot(potEnrgVector)
    end
    legend(cfg.j2spar.segmentName)
end
tsvFile.info.jointname = cfg.j2spar.segmentName;
tsvFile.processing.potEnerg = potEnrgVector;
end

