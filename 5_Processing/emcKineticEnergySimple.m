function [ tsvFile ] = emcKineticEnergySimple( tsvFile, cfg )
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
%     *.markerNum: vector contraining the markers index constituting all the segments
%     *.markerName: cell vector contraining the names of all the segments 
%     *.demsterIdx: vector contraining the index of the weights of all the
%     segments in the DEMSPTER table
%     [OPTIONAL]
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     *.featMarker: markers on which the feature is calculated
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.markerNum = {[4,9,11,12],[4,11],[4,5,1],[1,2,3]};
% cfg.markerName = {'PELVIS', 'RPELV', 'RTHIGH', 'RTIB'};
% cfg.demsterIdx = [0 1 8 7];
% cfg.paramComputation = 'mean';
% cfg.display = true;
% cfg.displayUnit = 's'
% tsvFile = emcKineticEnergySimple(tsvFile, cfg);
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

errorIfNotField(cfg, 'markerNum') % Transform markers to joints
errorIfNotField(cfg, 'markerName') % Transform joints to segments
errorIfNotField(cfg, 'demsterIdx') % Weights of segments



%% COMPUTATION AREA
disp('[PROCESSING] Kinetic Energy')
masssubject = 70; % Mass in kg
% Compute Speed
cfg.euclidianFlag = true;
cfg.deriv = 1;
tsvFileSpeed = emcDerivative(tsvFile, cfg);
speedMatrix = tsvFileSpeed.processing.speed;

% Compute mass
dempsterParam = mcgetsegmpar('Dempster', cfg.demsterIdx);
mass = dempsterParam.m*masssubject;

% Compute Kinetic Energy
kinEnerg = zeros(size(speedMatrix,1), numel(cfg.markerName));
for i = 1:numel(cfg.markerName)
    segmentName = cfg.markerName{i};
    segmentNum = cfg.markerNum{i};
    segmentMass = mass(i);
    segmentSpeed = speedMatrix(:,segmentNum);
    segmentKinEnergAll = 0.5*segmentMass*segmentSpeed.^2;
    kinEnerg(:,i) = mean(segmentKinEnergAll,2);
end

% Print the continuous values
if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Kinetic Energy'],'NumberTitle','off')
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        lineStyles = linspecer(size(kinEnerg,2));
%         lineStyles = distinguishable_colors(size(kinEnerg,2));
        hold on
        for i = 1:size(kinEnerg,2)
            plot(time, kinEnerg(:,i),'Color',lineStyles(i,:))
        end
        hold off
    else
        lineStyles = linspecer(size(kinEnerg,2));
%         lineStyles = distinguishable_colors(size(kinEnerg,2));
        hold on
        for i = 1:size(kinEnerg,2)
            plot(kinEnerg(:,i),'Color',lineStyles(i,:))
        end
        hold off
    end
    legend(cfg.markerName)
end
tsvFile.info.segmentName = cfg.markerName;
tsvFile.processing.kinEnerg = kinEnerg;
end

