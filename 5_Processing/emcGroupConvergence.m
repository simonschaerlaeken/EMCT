function [ varargout ] = emcGroupConvergence( varargin )
% Computes convergence of lines of sight. Either globally between a fixed
% point and the intersection of all lines for sight or partially between the
% different bodies
% 
% syntax
% tsvFile = emcGroupConvergence(tsvFile,tsvFile2[, tsvFile3,...], cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.convergenceType: str type of convergence - global: calculating the
%     distance between the intersection of all line of sight and a fixed
%     point. Partial: calculating the distance between one body and the
%     line of sight converging of all the other bodies
%     *.headMarker: cell array containing the markernames for the head
%     defininf the line of sight (suggestion: use baricenter and middlepoint)
%     [OPTIONAL]
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     *.displayVerification: boolean deciding if a figure is to be plotted  for verification of the line of sights(default: false)
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.convergenceType = 'global';
% cfg.headMarker = {'HEADBar','HEADMid'}
% cfg.displayUnit = 's'
% cfg.display = true
% tsvFile = emcGroupConvergence(tsvFile,tsvFile2, tsvFile3, cfg));
% 
% comments
% feature is saved in tsvFile.processing.groupConvGlobal/tsvFile.processing.groupConvPartial
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
if numel(varargin) < 2
    error('Not enough argument, Min: tsvFile1, tsvFile2, cfg')
end
cfg = varargin{end}; % CFG is the last argument
tsvFileIdxList = 1:numel(varargin)-1;
varargout = varargin(1:end-1);

for tsvFileIdx = 1:numel(tsvFileIdxList)
    if ~isfield(varargin{tsvFileIdx}, 'info')||~isfield(varargin{tsvFileIdx}.info, 'filename')
        filename = strsplit(varargin{tsvFileIdx}.filename, filesep);
        filename = filename{end};
        varargin{tsvFileIdx}.info.filename = filename;
    end
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
% display - if true, plot the verification of the intersection of the line
% of sights
if ~isfield(cfg, 'displayVerification')
    disp('cfg.displayVerification was not setup. Default: false')
    cfg.displayVerification = false;
end


% convergenceType - global or partial
errorIfNotField(cfg, 'convergenceType')

% Check error
errorIfNotField(cfg, 'headMarker') % Marker for the head



%% COMPUTATION AREA
% Intitialise the vector containing the 2 points for the head for all the
% bodies (2D)
points = zeros(size(varargin{1}.data,1),numel(tsvFileIdxList)*4); % *4 because you don't need the Z coordinate
% Retrieve all the points for the vector for the head (for each body)
% Use only 2D information (get rid of Z-values)
for tsvFileIdx = 1:numel(tsvFileIdxList)
    tsvFile = varargin{tsvFileIdx};
    % Extract data for points of interest
    headMarkerIdxList = findIndexList(tsvFile.markerName, cfg.headMarker);
    tsvFileMarker = mcgetmarker(tsvFile, headMarkerIdxList);
    dataPoints = tsvFileMarker.data(:,[1 2 4 5]); % Delete z coordinate
    % Insert in points
    pointIdx = ((tsvFileIdx-1)*4)+1;
    points(:,pointIdx:pointIdx+3) = dataPoints;
end
% Define the permutation to check the intersection of two lines for every
% possibilities
permutationBody = nchoosek(1:numel(tsvFileIdxList),2);
% Initialise the intersection points
intersectionPoints = zeros(size(points,1),size(permutationBody,1)*2);
% Compute one permutation at the time
for i = 1:size(permutationBody,1)
    % Select two bodies
    bodyIdx1 = permutationBody(i,1);
    bodyIdx2 = permutationBody(i,2);
%     disp(sprintf('Permutation = %i,%i', body_idx_1,body_idx_2));
    % Extract their respective points
    pointBodyIdx = ((bodyIdx1-1)*4)+1;
    pointsBody1 = points(:,pointBodyIdx:pointBodyIdx+3);
    pointBodyIdx = ((bodyIdx2-1)*4)+1;
    pointsBody2 = points(:,pointBodyIdx:pointBodyIdx+3);
    % Find the intersection between the two lines for every frame
    intersectionIdx = ((i-1)*2)+1;
    [intersectionTmp, ~] = findIntersection2D(pointsBody1, pointsBody2);
    intersectionPoints(:,intersectionIdx:intersectionIdx+1) = intersectionTmp;
end
% PARTIAL / GLOBAL
if strcmp('global',cfg.convergenceType)
    % Find the ear point (reference point to mesure distance with baricentre)
    if ~isfield(cfg, 'earPoint') || isempty(cfg.earPoint)
        cfg.earPoint = getEarPointGraph(points);
        cfg.earPoint = repmat(cfg.earPoint, size(points,1),1);
    else
        % Previously calculated ear_point does not have the same size
        if size(cfg.earPoint,1) ~= size(points,1)
            cfg.earPoint = repmat(cfg.earPoint(1,:), size(points,1),1);
        end
    end
    baricentrePoint = baricentre2D(intersectionPoints);
    if cfg.display && cfg.displayVerification
        title = 'Global Convergence - Baricenter, intersections, and ear point';
        verificationPlot(cfg, points, intersectionPoints, baricentrePoint, cfg.earPoint, title)
    end
    % Compute the distance between the point of conv and the ear point
    convGlobal = zeros(size(baricentrePoint,1),1);
    for i = 1:size(baricentrePoint,1)
        convGlobal(i) = pdist([baricentrePoint(i,:);cfg.earPoint(i,:)], 'euclidean');
    end
    if cfg.display
        figure('Name', 'Distance between a fixed point and baricentre of intersections','NumberTitle','off')
        if strcmp(cfg.displayUnit, 's')
            time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
            plot(time, convGlobal)
        else
            plot(convGlobal)
        end
        %trimmean_dist = trimmean(distance_bar_ear,cfg_main.group_trimmean_param);
        hold on
        trmDist = trimmean(convGlobal,30);
        plot (xlim, [trmDist trmDist], 'r');
        text(0,trmDist,'\leftarrow trimmean(dist)')
        hold off
        legend('Euclidian Distance')
    end
    % Modify the tsvFiles
    for tsvFileIdx = 1:numel(tsvFileIdxList)
        varargout{tsvFileIdx}.processing.groupConvGlobal = convGlobal;
    end
% PARTIAL
elseif strcmp('partial',cfg.convergenceType)
    convPartial = [];
    for tsvFileIdx = 1:numel(tsvFileIdxList)
        % Find baricenter for body of interest
        pointIdx = ((tsvFileIdx-1)*4)+1;
        pointsBody = points(:,pointIdx:pointIdx+3);
        % Find baricenter for the intersection of the line of sight of the
        % other bodies
        indexSubsetPermutation = findSubsetPermutation(permutationBody,i);
        subsetIntersection = intersectionPoints(:,indexSubsetPermutation);
        baricentrePartialPoint = baricentre2D(subsetIntersection);
        % Display
        if cfg.display && cfg.displayVerification
            title = ['Partial Convergence - Baricenter, intersections, body of interest =' tsvFileIdx];
            verificationPlot(cfg, points, subsetIntersection, baricentrePartialPoint, pointsBody(:,1:2), title)
        end
        % Compute the distance between the point of conv and the ear point
        distancePartial = zeros(size(baricentrePartialPoint,1),1);
        for t = 1:size(baricentrePartialPoint,1)
            distancePartial(t) = pdist([baricentrePartialPoint(t,:);pointsBody(t,1:2)], 'euclidean');
        end
        convPartial = horzcat(convPartial,distancePartial);
    end
    if cfg.display
        figure('Name', 'Distance between a partial point of convergence and respective left-out body','NumberTitle','off')
        if strcmp(cfg.displayUnit, 's')
            time = 0:(1/tsvFileMarker.freq):((tsvFileMarker.nFrames-1)/tsvFileMarker.freq);
            plot(time, convPartial)
        else
            plot(convPartial)
        end
        % Legend
        % Modify each tsvFile
        legendValues = cell(1, numel(tsvFileIdxList));
        for tsvFileIdx = 1:numel(tsvFileIdxList)
             legendValues{tsvFileIdx} = varargout{tsvFileIdx}.info.filename;
        end
        legend(legendValues)
    end
    % Modify each tsvFile
    for tsvFileIdx = 1:numel(tsvFileIdxList)
        varargout{tsvFileIdx}.processing.groupConvPartial = convPartial;
    end
end
end


% >Function to find the index of the permutation that dont contain a
% specific body
function [indexExpanded]= findSubsetPermutation(permutation,bodyIdx)
index = [];
for i = 1:size(permutation,1)
    if permutation(i,1) ~= bodyIdx && permutation(i,2) ~= bodyIdx
        index(end+1) = i;
    end
end
% expand
indexExpanded = [];
for i = 1:length(index)
    indexTmp = ((index(i)-1)*2)+1;
    indexExpanded(end+1:end+2) = [indexTmp, indexTmp+1];
end
end

% Function to print the graph with the head position and choose a location
% for the ear point
function [earPoint]= getEarPointGraph(points)
% Generate random point in time
a = 0;
b = size(points,1);
randomNumber = int64(round((b-a).*rand(1,1) + a));
points = points(randomNumber,:);
% Reshape points 
pointsReshaped = reshape(points, [2,length(points)/2]);
pointsReshaped = pointsReshaped';
% Create Line
lines = [];
for i = 1:2:size(pointsReshaped,1)
    lines = vertcat(lines, createLine(pointsReshaped(i,:), pointsReshaped(i+1,:)));
end
% prepare a figure for drawing
% draw the triangle
figure; clf;
hold on; 
% axis([0 20 0 20]);
axis equal;
drawPoint(pointsReshaped, 'marker', 'o', 'markersize', 8, 'linewidth', 2, ...
    'markerFaceColor', 'w');
drawLine(lines, 'color', 'r');
earPoint = ginput(1);
close
% close clf
end


% Finds the point of intersection of two lines in 2D space.
% Lines are defined by two points.
function [P, flag]= findIntersection2D(points1, points2)
% nd = 2;
nf = size(points1, 1);
% Initialize
P = zeros(nf, 2);
t = zeros(nf, 2);
flag = zeros(nf, 1);
% Compute for every frame
for f = 1:nf
    [P(f,:), t(f,:)] = lineintersect(points1(f,:), points2(f,:));
    if (isnan(P(f,1)) || t(f,1) <= 0 || t(f,2) <= 0)
       flag(f,1) = 0;
    else
        flag(f,1) = 1;
    end
end
end


% Finds 2D intersection of two lines.
% Points and vectors are to be defined as row vectors.
% Returns NaN if the lines are parallel.
function [P, t] = lineintersect(points1, points2)
% Create vectors from 2 points
v1 = [(points1(3)-points1(1)) (points1(4)-points1(2))];
v2 = [(points2(3)-points2(1)) (points2(4)-points2(2))];
nd = size(points1)/2;
% Formula
b = points2(1:2)' - points1(1:2)';
A = [v1', -v2'];

if rank(A) < nd
    t = [NaN, NaN];
    P = [NaN, NaN];
else
    t = A\b;
    t = t';
    P = points1(1:2) + t(1)*v1;
end
end

% Finds 2D intersection of two lines.
% Points and vectors are to be defined as row vectors.
% Returns NaN if the lines are parallel.
function [bariPoints] = baricentre2D(points)
% Sum all the points
bariPoints(:,1) = nansum(points(:,1:2:size(points,2)),2);
bariPoints(:,2) = nansum(points(:,2:2:size(points,2)),2);
% Divide them by nb of points
bariPoints = bariPoints/(size(points,2)/2);
end

% Function to plot in 3D the representation of the point and lines
% ----------------------------------------------------------------
function [] = verificationPlot(cfg, points, intersectionPoints, baricentrePoint, earPoint, title)
% Frame
if ~isfield(cfg, 'frame') || isempty(cfg.frame)
    % Generate random point in time
    a = 1;
    b = size(points,1);
    frame = int64(round((b-a).*rand(1,1) + a));
else
    frame = cfg.frame;
end
% Create Points
points = points(frame,:);
% Reshape points 
pointsReshaped = reshape(points, [2,length(points)/2]);
pointsReshaped = pointsReshaped';
% Create Lines
lines = [];
for i = 1:2:size(pointsReshaped,1)
    lines = vertcat(lines, createLine(pointsReshaped(i,:), pointsReshaped(i+1,:)));
end

% prepare a figure for drawing
% draw the triangle
figure('Name', ['At timeframe = ' num2str(frame) ', ' title],'NumberTitle','off'); clf;
hold on; 
% axis([0 20 0 20]);
axis equal;
drawPoint(pointsReshaped, 'marker', 'o', 'markersize', 8, 'linewidth', 2, ...
    'markerFaceColor', 'w');
drawLine(lines, 'color', 'r');

% Extra point
if ~isempty(intersectionPoints) && ~isempty(baricentrePoint) && ~isempty(earPoint)
    intersectionPoints = intersectionPoints(frame,:);
    baricentrePoint = baricentrePoint(frame,:);
    earPoint = earPoint(frame,:);

    % Reshape intersection points
    intersectionPointsReshaped = reshape(intersectionPoints, [2,length(intersectionPoints)/2]);
    intersectionPointsReshaped = intersectionPointsReshaped';

    drawPoint(intersectionPointsReshaped, 'marker', 'x', 'markersize', 8, 'linewidth', 2, ...
        'markerFaceColor', 'b');
    drawPoint(baricentrePoint, 'marker', 'x', 'markersize', 10, 'linewidth', 2, ...
        'color', 'g');
    drawPoint(earPoint, 'marker', 'o', 'markersize', 10, 'linewidth', 2, ...
        'color', 'g');
end
hold off;
end


