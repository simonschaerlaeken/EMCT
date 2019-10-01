function [ tsvFile ] = emcAngleVector( tsvFile, cfg )
% Computes the angle between 2 lines, each of which is defined by 2 markers
% 
% syntax
% tsvFile = emcAngleVector(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.angleName: str name of the angle
%     *.angleVector: cell array defining the 2 lines (4 markers)
%     [OPTIONAL]
%     *.angleDim: int value representing the dimension in which the angle
%     is calculated -- if 2: angle is computed only in the X-Y space. if 3:
%     angle is computed in the X-Y-Z space (default: 3)
%     feature is calculated (default: all)
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.angleName = 'bending';
% cfg.angleVector = {{'a','b'},{'c','d'}};
% cfg.angleDim = 3;
% cfg.display = true;
% cfg.displayUnit = 's'
% tsvFile = emcAngleVector(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.(angleName)
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

% angleName - define the name for the variable
errorIfNotField(cfg,'angleName')

% angleVector - define the 2 vectors composing the angle
errorIfNotField(cfg,'angleVector')
if numel(cfg.angleVector) ~= 2
    error('cfg.angleVector must be represented by a list of 2 cells containing the markers for each vector. Ex {{"a","b"},{"b","c"}}')
elseif numel(cfg.angleVector{1}) ~= 2 || numel(cfg.angleVector{1}) ~= 2
    error('cfg.angleVector must be represented by a list of 2 cells containing the 2 markers for each vector. Ex {{"a","b"},{"b","c"}}')
end
% angleDim - define how many dimension should be use to calculate the
% angle, can either be 2 or 3
if ~isfield(cfg, 'angleDim')
    disp('cfg.angleDim was not setup. Default: 3D')
    cfg.angleDim = 3;
elseif cfg.angleDim ~= 2 && cfg.angleDim ~= 3
    error('cfg.angleDim can only take value equals to 2 or 3')
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
%% COMPUTING AREA
% Vector 1
vector1MarkerIdx = findIndexList(tsvFile.markerName, cfg.angleVector{1});
tsvFileVector1 = mcgetmarker(tsvFile,vector1MarkerIdx);
% Vector 2
vector2MarkerIdx = findIndexList(tsvFile.markerName, cfg.angleVector{2});
tsvFileVector2 = mcgetmarker(tsvFile,vector2MarkerIdx);
% Concat
vectorTable = horzcat(tsvFileVector1.data, tsvFileVector2.data);
% Compute angle
angleTable = anglesWith2vectors(cfg, vectorTable);

% Display
if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Angle - ', cfg.angleName],'NumberTitle','off')
    subplot(1,2,1)
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        plot(time, angleTable)
    else
        plot(angleTable)
    end
    subplot(1,2,2)
    circ_plot(circ_ang2rad(angleTable),'hist',[],20,true,true,'linewidth',2,'color','r');
end

tsvFile.processing.(cfg.angleName) = angleTable;

end

%FEAT_ANGLES Angles between different vectors within the human body
function [angleTable] = anglesWith2vectors(cfg, vectorTable)
% Create the two vector from the data
v1 = vectorTable(:,4:6) - vectorTable(:,1:3);
v2 = vectorTable(:,10:12) - vectorTable(:,7:9);
% Project the vectors if only two dimensions
if cfg.angleDim == 2
    v1 = v1(:,1:2);
    v2 = v2(:,1:2);
    % Calcul Angle
    angleTable = zeros(size(vectorTable,1),1);
    for t = 1:size(vectorTable,1)
        % angle_table(t) = rad2deg(atan2(v2(t,2), v2(t,1)) - atan2(v1(t,2), v1(t,2)));
%         angle_table(t) = rad2deg(atan2(v2(t,2)- v1(t,2), v2(t,1)- v1(t,1)));
        angleTable(t) = acosd(dot(v1(t,:),v2(t,:))/(norm(v1(t,:))*norm(v2(t,:))));
        if (angleTable(t) > 90)
            angleTable(t) = 180 - angleTable(t);
        end
    end
elseif cfg.angleDim == 3
    % Calcul Angle
    angleTable = zeros(size(vectorTable,1),1);
    for t = 1:size(vectorTable,1)
        angleTable(t) = rad2deg(atan2(norm(cross(v1(t,:),v2(t,:))),dot(v1(t,:),v2(t,:))));
    end
end


% Circular Statistic to extract mean and variance
% angle_table_final = zeros(1,2);
% angle_table_rad = circ_ang2rad(angle_table);
% angle_table_rad_mean = circ_mean(angle_table_rad);
% angle_table_final(1) = rad2deg(angle_table_rad_mean);
% angle_table_rad_var = circ_var(angle_table_rad);
% angle_table_final(2) = rad2deg(angle_table_rad_var);
% angle_table = angle_table_final;

end

