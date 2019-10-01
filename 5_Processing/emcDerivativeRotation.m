function [ tsvFile ] = emcDerivativeRotation(tsvFile, cfg)
% Computes the derivative of the motion 
% 
% syntax
% tsvFile = emcDerivative(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.deriv: int order of the derivative to be calculated - 1 (speed), 2
%     (acceleration), 3 (jerk)
%     *.euclidianFlag: boolean to determine if the script uses euclidian
%     norm after derivative
%     [OPTIONAL]
%     *.featMarker: cell array containing the markernames on which the
%     feature is calculated (default: all)
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.featMarker = {'a','b','c'};
% cfg.deriv = 1;
% cfg.euclidianFlag = true;
% cfg.display = true
% cfg.displayUnit = 's'
% tsvFile = emcDerivative(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.speed
% feature is saved in tsvFile.processing.acceleration
% feature is saved in tsvFile.processing.jerk
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland

%% CHECKING AREA
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
errorIfNotField(tsvFile, 'rotationdata') % If euclidian norm or raw data as output

% featMarker - markers on which the feature is calculated
if ~isfield(cfg, 'featMarker')||isempty(cfg.featMarker) 
    disp('cfg.featMarker was not setup. Default: all markers')
    cfg.featMarker = tsvFile.rotationMarkerName;
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
% Check error
errorIfNotField(cfg, 'deriv') % Order of derivative
errorIfNotField(cfg, 'euclidianFlag') % If euclidian norm or raw data as output


%% COMPUTATION AREA
disp('[PROCESSING] Derivative')
% Transform Rotation into Position (So that we can use all the preestablished functions)
tsvFileRot = tsvFile;
rotationMarkerName = {};
for i = 1:numel(tsvFile.rotationMarkerName)
    for j = 1:tsvFile.rotationDim
        rotationMarkerName{end+1} = [tsvFile.rotationMarkerName{i}, '_', num2str(j)];
    end
end
tsvFileRot.markerName = rotationMarkerName;
tsvFileRot.nMarkers = numel(tsvFileRot.markerName);
tsvFileRot.data = tsvFile.rotationdata;
% % Select markers
% featMarkerIndex = findIndexList(tsvFile.rotationMarkerName, cfg.featMarker);
% featMarkerIndexRot = [];
% for i = 1:numel(featMarkerIndex)
%     for j = 0:2
%         featMarkerIndexRot(end+1) = ((featMarkerIndex(i)-1)*3+1)+(j*3);
%     end
% end
% tsvFileMarker = mcgetmarker(tsvFileRot, featMarkerIndexRot);
tsvFileMarker = tsvFileRot;
% Derivation
warning('off','all')
tsvFileDer = mctimeder(tsvFileMarker, cfg.deriv);
warning('on','all')
% Norm
tsvFileDerNorm = mcnorm(tsvFileDer);
derivMatrixNorm = tsvFileDerNorm.data;
if cfg.euclidianFlag
    tsvFileDerNorm.type = 'MoCap data';
    tsvFileDerNormNorm = mcnorm(tsvFileDerNorm);
    derivMatrixNorm = tsvFileDerNormNorm.data;
end

%% SMOOTHING
for i = 1:size(derivMatrixNorm,2)
    derivMatrixNorm(:,i) = smooth(derivMatrixNorm(:,i), 100);
end
%%
% Right title
if cfg.deriv == 1
    derivName = 'speed';
elseif cfg.deriv == 2
    derivName = 'acceleration';
elseif cfg.deriv == 3
    derivName = 'jerk';
end
% Display
if cfg.display
    colorset = distinguishable_colors(length(cfg.featMarker));
    % Add it to the title and then figure
    title = [tsvFile.info.filename, ' - ', derivName];
    figure('Name',title,'NumberTitle','off')
    if cfg.euclidianFlag
        hold on
        for i = 1:size(derivMatrixNorm,2)
            if strcmp(cfg.displayUnit, 's')
                time = 0:(1/tsvFileMarker.freq):((tsvFileMarker.nFrames-1)/tsvFileMarker.freq);
                plot(time, derivMatrixNorm(:,i), 'color', colorset(i,:))
            else
                plot(derivMatrixNorm(:,i), 'color', colorset(i,:))
            end
        end
        hold off
        legend(cfg.featMarker)
    else
        for i = 1:3:size(derivMatrixNorm,2)
            if strcmp(cfg.displayUnit, 's')
                time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
                hold on,ax1 = subplot(3,1,1);plot(time,derivMatrixNorm(:,i), 'color', colorset(i,:)),title('x');
                hold on,ax2 = subplot(3,1,2);plot(time,derivMatrixNorm(:,i+1), 'color', colorset(i,:)),title('y');
                hold on,ax3 = subplot(3,1,3);plot(time,derivMatrixNorm(:,i+2), 'color', colorset(i,:)),title('z');
            else
                hold on,ax1 = subplot(3,1,1);plot(derivMatrixNorm(:,i), 'color', colorset(i,:)),title('x');
                hold on,ax2 = subplot(3,1,2);plot(derivMatrixNorm(:,i+1), 'color', colorset(i,:)),title('y');
                hold on,ax3 = subplot(3,1,3);plot(derivMatrixNorm(:,i+2), 'color', colorset(i,:)),title('z');
            end
        end
        hold off
        legend(cfg.featMarker)
    end
end
% Define the correct output
if cfg.euclidianFlag
    derOutput = derivMatrixNorm;
else
    derOutput = derivMatrixNorm;
end
% Set up the output to the right structure
eval(['tsvFile.processing.',derivName,' = derOutput;'])

end


