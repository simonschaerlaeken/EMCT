function [ tsvFile ] = emcKineticEnergy( tsvFile, cfg )
% Computes the kinetic energy of the markers (based on mcpotenergy from the
% Motion Capture Toolbox, Copyright 2008, University of Jyvaskyla, Finland
% )
% 
% syntax
% tsvFile = emcKineticEnergy(tsvFile, cfg)
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
% tsvFile = emcKineticEnergy(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.kinEnerg
% Segment number values for model 'Dempster' are as follows: no parameter=0, 
% hand=1, forearm=2, upper arm=3, forearm and hand=4, upper extremity=5, 
% foot=6, leg=7, thigh=8, lower extremity=9, head=10, shoulder=11, thorax=12, 
% abdomen=13, pelvis=14, thorax and abdomen=15, abdomen and pelvis=16, 
% trunk=17, head, arms and trunk (to glenohumeral joint)=18, 
% head, arms and trunk (to mid-rib)=19
% see also
% mckinenergy
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



%% COMPUTATION AREA
disp('[PROCESSING] Kinetic Energy')
% Create Joint and segment files
tsvFileJoint = mcm2j(tsvFile, cfg.m2jpar);
tsvFileJoint.data = tsvFileJoint.data*1000; % If you save it it is too small -- don't understand why - maybe consider in mm while we are in meters
% emcSaveTSV2C3D(tsvFileJoint,[pwd,filesep, 'Test.c3d'])
%cfg.connectionMatrix = [1,5;1,2;1,3;5,4;4,6;4,7;7,8;8,9;9,10;11,4;11,12;13,14;15,2;15,16;16,17;18,3;18,19;19,20];
%emcPlotBody3D(tsvFileJoint, cfg);
tsvFileSegment = mcj2s(tsvFileJoint, cfg.j2spar);

% Compute Kinetic energy
kinEnrgVector = mckinenergy(tsvFileJoint, tsvFileSegment, cfg.mcspar);

% Print the continuous values
if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Kinetic Energy'],'NumberTitle','off')
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        lineStyles = linspecer(size(kinEnrgVector,2));
        lineStyles = distinguishable_colors(size(kinEnrgVector,2));
        hold on
        for i = 1:size(kinEnrgVector,2)
            plot(time, kinEnrgVector(:,i),'Color',lineStyles(i,:))
        end
        hold off
    else
        lineStyles = linspecer(size(kinEnrgVector,2));
        lineStyles = distinguishable_colors(size(kinEnrgVector,2));
        hold on
        for i = 1:size(kinEnrgVector,2)
            plot(kinEnrgVector(:,i),'Color',lineStyles(i,:))
        end
        hold off
    end
    legend(cfg.j2spar.segmentName)
end
tsvFile.info.jointname = cfg.j2spar.segmentName;
tsvFile.processing.kinEnerg = kinEnrgVector;
end

