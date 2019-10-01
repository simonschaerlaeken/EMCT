function [] = emcPlotCompare(varargin)
% Plots the motion capture data of markers selected on the
% three axes X, Y, and Z. Can take as many tsvfile as wanted and plot them
% in different colors to compare them.
% 
% syntax
% emcPlotCompare(tsvFile1[, tsvFile2,...], cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.compareType: string type of comparison. 'marker', 'feature'
%     IF compareType: 'feature'
%       *.featureList: cell array containing the feature names to be plotted
%       *.featureAverageFlah: boolean to specify if you want the average of
%     the feature to appear on the plot
%     [OPTIONAL]
%     *.markerList: cell array with markers to be plotted (if not defined all markers will be used)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
%     
% output
% -
% 
% examples
% cfg.compareType = 'marker'
% cfg.markerList = {'a','b','c'};
% cfg.displayUnit = 's'
% emcPlotCompare(tsvFile1, tsvFile2, cfg);
%
% cfg.compareType = 'feature'
% cfg.featureList = {'speed'};
% cfg.displayUnit = 's'
% emcPlotCompare(tsvFile1, tsvFile2, cfg);
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
if nargin<2
    error('ERROR: No file or no cfg as input')
elseif nargin<3
    warning('Only one file to compare, last argument is cfg')
end
% cfg in last argument
cfg = varargin{end};

% compareType - string type of comparison. 'marker', 'feature'
errorIfNotField(cfg, 'compareType');
% Check for each arg
for fileIdx = 1:nargin-1
    % filename - for title
    if ~isfield(varargin{fileIdx}, 'info')||~isfield(varargin{fileIdx}.info, 'filename')
        filename = strsplit(varargin{fileIdx}.filename, '\');
        filename = filename{end};
        varargin{fileIdx}.info.filename = filename;
    end
    if strcmp(cfg.compareType, 'marker')
        % plotMarkerList - choice of marker to display
        if ~isfield(cfg, 'plotMarkerList')||isempty(cfg.plotMarkerList)
            warning('No marker was define, by default all markers are printed. Define markers to be printed using "cfg.plotMarkerList"')
            cfg.plotMarkerList = varargin{fileIdx}.markerName;
        end
    elseif strcmp(cfg.compareType, 'feature')
        errorIfNotField(cfg, 'featureList');
        for featureIdx = 1:numel(cfg.featureList)
            featureName = cfg.featureList{featureIdx};
            if ~isfield(varargin{fileIdx}.processing, featureName)
                error(['Feature ', featureName, 'does not exist in tsvFile.processing'])
            end
        end
    else
        error('cfg.compareType can either be "marker" or "feature"');
    end
end
% featureAverageFlag: boolean to specify if you want the average of the feature to appear on the plot
if ~isfield(cfg, 'featureAverageFlag')
    cfg.featureAverageFlag = false;
end
% displayUnit - str containing the unit of duration - 's' or 'tf'
if ~isfield(cfg, 'displayUnit')
    disp('cfg.displayUnit was not setup. Default: "tf"')
    cfg.displayUnit = 'tf';
end
%% COMPUTATION AREA
% Name to compare:
filename = '';
filenameList = cell(nargin-1,1);
for fileIdx = 1:nargin-1
    filename = strcat(filename, varargin{fileIdx}.info.filename, '-');
    filenameList{fileIdx} = varargin{fileIdx}.info.filename;
end
disp(['[PLOT] ', filename ' - Compare TSV']);
if strcmp(cfg.compareType, 'marker')
    % Define marker_selected with prefix for body
    hFig = figure('Name',[filename ' -- Compare Markers'], 'NumberTitle','off'); clf; hold on;
    set(hFig, 'Position' , [230 230 1300 700])

    % Create a different color for each value
    colors_individual = distinguishable_colors(nargin);

    % Process the different bodies
    for fileIdx = 1:nargin-1
        % Extract data for specific marker
        tsvFile = varargin{fileIdx};
        markerSelectedIdx = findIndexList(tsvFile.markerName, cfg.plotMarkerList);
        tsvFileMarker = mcgetmarker(tsvFile, markerSelectedIdx);
        %time = 1:length(tsv_file_marker.data);
        for i = 1:3:length(tsvFileMarker.markerName)*3
            % Print Graph
            if strcmp(cfg.displayUnit, 's')
                time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
                hold on,ax1 = subplot(3,1,1);plot(time, tsvFileMarker.data(:,i),'color',colors_individual(fileIdx,:)),title('x')
                hold on,ax2 = subplot(3,1,2);plot(time, tsvFileMarker.data(:,i+1),'color',colors_individual(fileIdx,:)),title('y')
                hold on,ax3 = subplot(3,1,3);plot(time, tsvFileMarker.data(:,i+2),'color',colors_individual(fileIdx,:)),title('z')
            else
                hold on,ax1 = subplot(3,1,1);plot(tsvFileMarker.data(:,i),'color',colors_individual(fileIdx,:)),title('x')
                hold on,ax2 = subplot(3,1,2);plot(tsvFileMarker.data(:,i+1),'color',colors_individual(fileIdx,:)),title('y')
                hold on,ax3 = subplot(3,1,3);plot(tsvFileMarker.data(:,i+2),'color',colors_individual(fileIdx,:)),title('z')
            end
        end
    end
    % Set global Title
    % annotation(h,'textbox',[0.49 0.97 0.12 0.032],...
    %     'String',title_fig,...
    %     'FitBoxToText','on',...
    %     'BackgroundColor',[1 1 1]);
elseif strcmp(cfg.compareType, 'feature')
    for featureIdx = 1:numel(cfg.featureList)
        % Feature Name
        featureName = cfg.featureList{featureIdx};
        % Create a different color for each value
        colors_individual = distinguishable_colors(nargin);
        % Figure
        figure('Name',[filename, ' - ', featureName],'NumberTitle','off');
        hold on
         % Process the different bodies
        for fileIdx = 1:nargin-1
            % Extract data for specific marker
            tsvFile = varargin{fileIdx};
            % Feature Data
            featureData = tsvFile.processing.(featureName);
            % Plot
            if size(featureData,1) == 1
                disp(['Feature ', featureName, 'for tsvFile ', tsvFile{fileIdx}.info.filename, ' ='])
                disp(featureData)
            else
                % if average display
                if cfg.featureAverageFlag
                    for i = 1:size(featureData,2)
                        if strcmp(cfg.displayUnit, 's')
                            time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
                            p1 = plot(time, featureData(:,i),'Color',colors_individual(fileIdx,:));
                        else
                            p1 = plot(featureData(:,i),'Color',colors_individual(fileIdx,:));
                        end
                        p1.Color(4) = 0.75;
                    end
                    if strcmp(cfg.displayUnit, 's')
                        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
                        plot(time, mean(featureData,2),'LineWidth',3, 'Color',colors_individual(fileIdx,:));
                    else
                        plot(mean(featureData,2),'LineWidth',3, 'Color',colors_individual(fileIdx,:));
                    end
                % if not, special plot with legend
                else
                    if strcmp(cfg.displayUnit, 's')
                        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
                        eval(['h',num2str(fileIdx),' = plot(time,featureData,''color'',colors_individual(fileIdx,:));']);
                    else
                        eval(['h',num2str(fileIdx),' = plot(featureData,''color'',colors_individual(fileIdx,:));']);
                    end
                end
            end
        end
        % Set legend if not average display
        if ~cfg.featureAverageFlag
            % Set Legend
            hString = '';
            for fileIdx = 1:nargin-1
                hString = [hString, ' h', num2str(fileIdx),'(1)'];
            end
            eval(['legend([',hString,'], filenameList)']);
        end
        hold off
    end
end

