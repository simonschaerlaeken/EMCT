function [space] = emcHeatmapData( varargin )
% Computes the heatmap of the baricenter motion seen only on a 2D plane (x-y) from above. 
% 
% syntax
% tsvFile = emcHeatmap(dataMotion[, dataFeature], cfg)
% 
% input parameters
% dataMotion: table motion data (3 columns per marker, each line represent the position of the marker at time t)
% dataFeature: table feature data
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.reduceFlag: boolean flag meaning that the final resolution will be
%     reduced
%     if reduceFlag: true
%       *.reduceResolution: [mandatory] table with 2 values - the reduction
%       factor for X and Y.
%     *.display: boolean deciding if a figure is to be plotted (default: true)
%     
% output
% -
% 
% examples
% cfg.reduceFlag = true;
% cfg.reduceResolution = [2,3];
% emcHeatmapData(motionData, featureData, cfg);
% 
% comments
% -
% 
% see also
% emcHeatmap
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
if numel(varargin) > 3
    error('Too many input arguments. Max:data, feature, cfg')
elseif numel(varargin) < 2
    error('Too few input arguments. min: data, cfg')
end

cfg = varargin{end};

if numel(varargin) == 2
    data = varargin{1};
    feature = [];
elseif numel(varargin) == 3
    data = varargin{1};
    feature = varargin{2};
end

if ~isfield(cfg, 'reduceFlag')
    cfg.reduceFlag = false;
else
    errorIfNotField(cfg, 'reduceResolution')
end

if ~isfield(cfg, 'display')
    cfg.display = true;
end

%% COMPUTATION AREA
% Create baricenter for a specific body in a specific file
baricenterPoint = centroid3D(data);
% Keep only X and Y coordinates
baricenterPoint = baricenterPoint(:,1:2);
% Print the scatter plot
if cfg.display
    figure('Name','Displacement','NumberTitle','off')
    a = 25;
    c = linspace(1,10,length(baricenterPoint(:,1)));
    scatter(baricenterPoint(:,1),baricenterPoint(:,2),a,c,'filled') 
end

% Round baricenter to fit within a matrix
baricenterPoint = round(baricenterPoint);

% Feature modulated
if isempty(feature)
    % Create the Heatmap
    [space, baricenterModified] = hmCreate(baricenterPoint, 10, []);
    figure('Name','Space Heatmap','NumberTitle','off')
    imagesc(space)
    hold on
    a = 10;
    c = linspace(1,10,length(baricenterModified(:,1)));
    scatter(baricenterModified(:,2),baricenterModified(:,1),a,c,'filled') 
    hold off
    view([-90 90])
    title('Space Heatmap')
else
    % Create the Heatmap
    [space, baricenterModified] = hmCreate(baricenterPoint, 10, feature);
    figure('Name',['Space Heatmap modulated by ' cfg.heatmapFeature],'NumberTitle','off')
    imagesc(space)
    hold on
    a = 10;
    c = linspace(1,10,length(baricenterModified(:,1)));
    scatter(baricenterModified(:,2),baricenterModified(:,1),a,c,'filled') 
    hold off
    view([-90 90])
    title(['Space Heatmap modulated by ' cfg.heatmapFeature '[Attention Coordinates]'])
end
% Reduce space
if cfg.reduceFlag
    resolution.x = cfg.reduceResolution(1);
    resolution.y = cfg.reduceResolution(2);
    reducedSpace = zeros(resolution.x, resolution.y);
    reducedStep.x = round(size(space,1)/resolution.x);
    reducedStep.y = round(size(space,2)/resolution.y);
    for x = 1:resolution.x
        xSpace = ((x-1)*reducedStep.x)+1;
        % Check if doesn't go above limits
        xSpaceFrame = xSpace+reducedStep.x;
        if xSpaceFrame > size(space,1)
            xSpaceFrame = size(space,1);
        end
        for y = 1:resolution.y
            ySpace = ((y-1)*reducedStep.y)+1;
            % Check if doesn't go above limits
            ySpaceFrame = ySpace+reducedStep.y;
            if ySpaceFrame > size(space,2)
                ySpaceFrame = size(space,2);
            end
            % Sum over the area
            reducedSpace(x,y) = sum(sum(space(xSpace:xSpaceFrame,ySpace:ySpaceFrame)));
        end
    end
    
    if cfg.display
        figure('Name','Reduced Space Heatmap','NumberTitle','off')
        imagesc(reducedSpace)
        view([-90 90])
        title('Reduced Space Heatmap [Attention Coordinates]')
    end
end


end






