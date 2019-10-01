function [ tsvFile ] = emcGraphTheorySetup(tsvFile, cfg)
% Set up the graph for a graph theory analysis - compute the edges values
% and define pairs
% 
% syntax
% tsvFile = emcSampleEntropy(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.edgeFunction: str function to compute the edge values - "distance",
%     "correlation", "coherence"
%     if cfg.edgeFunction ~=  "distance"
%        *.featureName: str feature name on which the edgeFunction will be
%         computed
%     [OPTIONAL]
%     *.GTpairlist: cell array defining a list of pairs of 2 markers (default: ALL)
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.edgeFunction = 'correlation';
% cfg.featureName = 'kinEnerg';
% cfg.display = true;
% cfg.GTpairlist = {{'Marker1','Marker2'},{'Marker1','Marker3'}};
% tsvFile = emcGraphTheorySetup(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.GTsetup
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

% display - if true, plot the result
if ~isfield(cfg, 'display')
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
    errorIfNotField(cfg, 'connectionMatrix')
else cfg.display % If true must check for connection matrix
    errorIfNotField(cfg, 'connectionMatrix')
end
% frame - frame at which the body will be plotted
if ~isfield(cfg, 'frame')||isempty(cfg.frame) % Must choose random frame
    % Generate random point in time
    a = 1;
    b = size(tsvFile.data,1);
    cfg.frame = int64(round((b-a).*rand(1,1) + a));    
end
% edgeFunction
errorIfNotField(cfg,'edgeFunction')
% Parameter edgeFunction
if ~strcmp(cfg.edgeFunction, 'distance')
    errorIfNotField(cfg, 'featureName')
    disp('[WARNING] If not with all markers, watch out for correspondance')
end
if ~isfield(cfg, 'GTpairlist')||isempty(cfg.GTpairlist)
    disp('cfg.GTpairlist was not setup. Default: All connections. Computing the list...')
    cfg.GTpairlist = {};
    for markerIdx = 1:numel(tsvFile.markerName)
        for markerIdx2 = (markerIdx+1):numel(tsvFile.markerName)
            cfg.GTpairlist{end+1} = {tsvFile.markerName{markerIdx},tsvFile.markerName{markerIdx2}};
        end
    end
elseif ischar(cfg.GTpairlist)&&strcmp(cfg.GTpairlist,'connectionMatrix')
    errorIfNotField(cfg, 'connectionMatrix')
    cfg.GTpairlist = emcGraphTheoryConnectionMatrix(tsvFile, cfg.connectionMatrix);
end

%% COMPUTATION AREA
disp('[PROCESSING] Graph Theory Setup')
% Element necessary to compute graph
source = {};
target = {};
edge = [];
if strcmp(cfg.edgeFunction, 'distance')
    source = {};
    target = {};
    edge = [];
    for pairIdx = 1:numel(cfg.GTpairlist)
        featMarkerIndex = findIndexList(tsvFile.markerName, cfg.GTpairlist{pairIdx});
        tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
        % Calculate distance between markers
        distance = zeros(size(tsvFileMarker.data,1),1);
        %Compute distance for everu timeframe
        for i = 1:size(tsvFileMarker.data,1)
            distance(i) = distancePoints3d(tsvFileMarker.data(i,1:3),tsvFileMarker.data(i,4:6));
        end
        edge(end+1) = mean(distance);
        source{end+1} = cfg.GTpairlist{pairIdx}{1};
        target{end+1} = cfg.GTpairlist{pairIdx}{2};
    end
    tsvFile.processing.graph = graph(source,target,edge);
    
    LWidths = 5*tsvFile.processing.graph.Edges.Weight/max(tsvFile.processing.graph.Edges.Weight);
    %plot(tsvFile.processing.graph,'EdgeLabel',tsvFile.processing.graph.Edges.Weight,'LineWidth',LWidths)
    plot(tsvFile.processing.graph,'LineWidth',LWidths)
    
elseif strcmp(cfg.edgeFunction, 'correlation')
    for pairIdx = 1:numel(cfg.GTpairlist)
        featMarkerIndex = findIndexList(tsvFile.markerName, cfg.GTpairlist{pairIdx});
        %Compute correlation for the two vectors
        correlation = corrcoef(tsvFile.processing.(cfg.featureName)(:,featMarkerIndex(1)),tsvFile.processing.(cfg.featureName)(:,featMarkerIndex(2)));
        edge(end+1) = correlation(1,2);
        source{end+1} = cfg.GTpairlist{pairIdx}{1};
        target{end+1} = cfg.GTpairlist{pairIdx}{2};
    end
    tsvFile.processing.graph = graph(source,target,edge);
    
    % LINE WIDTH
    LWidths = 5*(tsvFile.processing.graph.Edges.Weight-min(tsvFile.processing.graph.Edges.Weight))/(max(tsvFile.processing.graph.Edges.Weight)-min(tsvFile.processing.graph.Edges.Weight))+1;
    % COLOR
    LColorsInd = 100*(tsvFile.processing.graph.Edges.Weight-min(tsvFile.processing.graph.Edges.Weight))/(max(tsvFile.processing.graph.Edges.Weight)-min(tsvFile.processing.graph.Edges.Weight));
    LColors = zeros(size(LColorsInd,1),3);
    cmap = colormap;
    for LColorsIndIdx = 1:length(LColorsInd)
        LColors(LColorsIndIdx,:) = ind2rgb(round(LColorsInd(LColorsIndIdx)),cmap);
    end
    % COORDINATE
    [XData, YData, ZData] = emcFindCoordinateGraph(tsvFile, cfg, tsvFile.processing.graph);
    % PLOT
    %plot(tsvFile.processing.graph,'EdgeColor',LColors,'EdgeLabel',tsvFile.processing.graph.Edges.Weight,'LineWidth',LWidths,'XData',YData,'YData',ZData)
    plot(tsvFile.processing.graph,'EdgeColor',LColors,'EdgeLabel',tsvFile.processing.graph.Edges.Weight,'LineWidth',LWidths,'XData',XData,'YData',YData,'ZData',ZData)
    %plot(tsvFile.processing.graph,'LineWidth',LWidths)
    
elseif strcmp(cfg.edgeFunction, 'coherence')
    warning('Work in progress')
end
end

function [GTpairlist] = emcGraphTheoryConnectionMatrix(tsvFile, connectionMatrix) 
GTpairlist = {};
for pairIdx = 1:size(connectionMatrix,1)
    marker1 = tsvFile.markerName{connectionMatrix(pairIdx,1)};
    marker2 = tsvFile.markerName{connectionMatrix(pairIdx,2)};
    GTpairlist{end+1} = {marker1,marker2};
end
end

function [XData, YData, ZData] = emcFindCoordinateGraph(tsvFile,cfg,G) 
listNodes = table2cell(G.Nodes);
featMarkerIndex = findIndexList(tsvFile.markerName, listNodes);
featMarkerIndexX = (featMarkerIndex-1)*3+1;
featMarkerIndexY = (featMarkerIndex-1)*3+2;
featMarkerIndexZ = (featMarkerIndex-1)*3+3;
XData = tsvFile.data(cfg.frame,featMarkerIndexX);
YData = tsvFile.data(cfg.frame,featMarkerIndexY);
ZData = tsvFile.data(cfg.frame,featMarkerIndexZ);
end