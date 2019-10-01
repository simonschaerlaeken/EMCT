function [ tsvFile ] = emcScramble( tsvFile, cfg )
% Spatially scramble the marker data 
% 
% syntax
% tsvFile = emcScramble( tsvFile, cfg );
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.scrambleType: string containing the type of scrambling
%     [OPTIONAL]
%     -
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.scrambleType = 'random';
% tsvFile = emcMidpointBaricenter(tsvFile, cfg);
% 
% comments
% The markers created are a combination of the markername and "Vert" for Vertical line
% This can be further use to compare angle to a straight line or to have an idea
% of how bended the movement is compared to a straight vertical line
% 
% see also
% emcMidpointBaricenter
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
% Check error
if ~isfield(cfg, 'scrambleType') % Marker to compute the baricenter and the midpoint
    cfg.scrambleType = 'random';
elseif strcmp(cfg.scrambleType, 'parts')
    errorIfNotField(cfg, 'scrambleParts')
end
%% COMPUTATION AREA
if strcmp(cfg.scrambleType,'random')
    % Change Name
    tsvFile.info.filename = [tsvFile.info.filename(1:end-4),'_RandScr.tsv'];
    
    minData = zeros(1,3);
    maxData = zeros(1,3);
    % Compute the extreme value for every axis
    xData = tsvFile.data(:,1:3:size(tsvFile.data, 2));
    minData(1) = min(min(xData));
    maxData(1) = max(max(xData));
    yData = tsvFile.data(:,2:3:size(tsvFile.data, 2));
    minData(2) = min(min(yData));
    maxData(2) = max(max(yData));
    zData = tsvFile.data(:,3:3:size(tsvFile.data, 2));
    minData(3) = min(min(zData));
    maxData(3) = max(max(zData));

    % Define a random starting point for each marker, compute the different for
    % each coordinate and apply it to the rest of the trajectory
    for markerIdx = 1:numel(tsvFile.markerName)
        %disp(tsvFile.markerName{markerIdx})
        markerIn = ((markerIdx-1)*3)+1;
        markerOut = markerIn + 2;
        % Marker traj
        minMarkerData = min(tsvFile.data(:,markerIn:markerOut));
        maxMarkerData = max(tsvFile.data(:,markerIn:markerOut));
        % Distance from start position
        startPos = tsvFile.data(1,markerIn:markerOut);
        distMinStart = abs(startPos - minMarkerData);
        distMaxStart = abs(maxMarkerData - startPos); 
        % Random position:
        %available = (maxData-distMaxStart)-(minData+distMinStart);
        available = (maxData)-(minData);
        %position = minData + distMinStart + (rand()*available);
        position = minData + (rand(1,3).*available);
        % Translation vector
        translationVect = position - tsvFile.data(1,markerIn:markerOut);
        for i = 1:size(tsvFile.data,1)
            tsvFile.data(i,markerIn:markerOut) = tsvFile.data(i,markerIn:markerOut)+translationVect;
        end
    end
end
if strcmp(cfg.scrambleType,'parts')
    % Change Name
    tsvFile.info.filename = [tsvFile.info.filename(1:end-4),'_Scr.tsv'];
    % For each part
    for partIdx = 1:numel(cfg.scrambleParts)
        listMarker = findIndexList(tsvFile.markerName, cfg.scrambleParts{partIdx});
        tsvFileMarker = mcgetmarker(tsvFile, listMarker);
        % Init
        minData = zeros(1,3);
        maxData = zeros(1,3);
        % Compute the extreme value for every axis
        xData = tsvFileMarker.data(:,1:3:size(tsvFileMarker.data, 2));
        minData(1) = min(min(xData));
        maxData(1) = max(max(xData));
        yData = tsvFileMarker.data(:,2:3:size(tsvFileMarker.data, 2));
        minData(2) = min(min(yData));
        maxData(2) = max(max(yData));
        zData = tsvFileMarker.data(:,3:3:size(tsvFileMarker.data, 2));
        minData(3) = min(min(zData));
        maxData(3) = max(max(zData));
        % Define a random starting point for each marker, compute the different for
        % each coordinate and apply it to the rest of the trajectory
        for markerIdx = 1:numel(listMarker)
            markerIdxSel = listMarker(markerIdx);
            %disp(tsvFile.markerName{markerIdx})
            markerSelIn = ((markerIdxSel-1)*3)+1;
            markerSelOut = markerSelIn + 2;
            % Marker traj
            minMarkerData = min(tsvFile.data(:,markerSelIn:markerSelOut));
            maxMarkerData = max(tsvFile.data(:,markerSelIn:markerSelOut));
            % Distance from start position
            startPos = tsvFile.data(1,markerSelIn:markerSelOut);
            distMinStart = abs(startPos - minMarkerData);
            distMaxStart = abs(maxMarkerData - startPos); 
            % Random position:
            %available = (maxData-distMaxStart)-(minData+distMinStart);
            available = (maxData)-(minData);
            %position = minData + distMinStart + (rand()*available);
            position = minData + (rand(1,3).*available);
            % Translation vector
            translationVect = position - tsvFile.data(1,markerSelIn:markerSelOut);
            for i = 1:size(tsvFile.data,1)
                tsvFile.data(i,markerSelIn:markerSelOut) = tsvFile.data(i,markerSelIn:markerSelOut)+translationVect;
            end
        end
    end
end
if strcmp(cfg.scrambleType,'symmetry')
    % Change Name
    tsvFile.info.filename = [tsvFile.info.filename(1:end-4),'_Sym.tsv'];
    % For each part
    for symIdx = 1:numel(cfg.scrambleSymmetry)
        listMarker = findIndexList(tsvFile.markerName, cfg.scrambleSymmetry{symIdx});
        tsvFileMarker = mcgetmarker(tsvFile, listMarker);
        % Init Deplacement
        deplacementMatrix = tsvFileMarker.data(2:end,:) - tsvFileMarker.data(1:end-1,:);
        % Inverse (Symmetry)
        deplacementMatrix = -deplacementMatrix;
        % Start from the initial point and add deplacements.
        for markerIdx = 1:numel(listMarker)
            markerIdxSel = listMarker(markerIdx);
            %disp(tsvFile.markerName{markerIdx})
            markerSelIn = ((markerIdxSel-1)*3)+1;
            markerSelOut = markerSelIn + 2;
            % Marker in out for deplacement
            markerIn = ((markerIdx-1)*3)+1;
            markerOut = markerIn + 2;
            % Add deplacement
            for i = 1:(size(tsvFile.data,1)-1)
                tsvFile.data(i+1,markerSelIn:markerSelOut) = tsvFile.data(i,markerSelIn:markerSelOut) + deplacementMatrix(i,markerIn:markerOut);
            end
        end
    end
end

