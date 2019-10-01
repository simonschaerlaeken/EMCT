function [ syncQ, syncq ] = emcEventSyncOld( varargin )
% Computes the synchronicity matrix for different markers, bodies or
% features as defined by Albordo, 2016
% 
% syntax
% tsvFile = emcEventSync(tsvFile1[, tsvFile2], cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.syncType: str type of synchrony - 'intra' (for one body, compares inside one body), 'inter'
%     (for multiple body, compares bodies one against the other)
%     *.dataType: str type of data - 'marker' (on the marker motion),
%     'feature' (on the computed features)
%     [OPTIONAL]
%     *.dataList: cell array containing the markernames or feature names on which the
%     sync is calculated (default: all)
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     
% output
% syncQ: table overall degree of sync
% syncq: table overall degree of delay
% 
% examples
% cfg.syncType = 'intra';
% cfg.dataType = 'marker';
% cfg.dataList = {'a','b','c'};
% cfg.display = true
% tsvFile = emcEventSync(tsvFile, cfg);
% 
% cfg.syncType = 'inter';
% cfg.dataType = 'feature';
% cfg.dataList = {'speed', 'convexhull'};
% cfg.display = true
% tsvFile = emcEventSync(tsvFile1, tsvFile2, cfg);
% 
% comments
% feature is saved in tsvFile.processing.directness
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
if numel(varargin) < 2,
    error('Not enough argument, Min: tsvFile, cfg')
end
% CFG
cfg = varargin{end}; % CFG is the last argument
tsvFileIdxList = 1:numel(varargin)-1;
tsvFileAll = varargin(1:end-1);
varargout = tsvFileAll;

for tsvFileIdx = 1:numel(tsvFileIdxList),
    if ~isfield(tsvFileAll{tsvFileIdx}, 'info')||~isfield(tsvFileAll{tsvFileIdx}.info, 'filename'),
        filename = strsplit(tsvFileAll{tsvFileIdx}.filename, filesep);
        filename = filename{end};
        tsvFileAll{tsvFileIdx}.info.filename = filename;
    end
end

% display - if true, plot the result
if ~isfield(cfg, 'display'),
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
end


% syncType - intra or inter
errorIfNotField(cfg, 'syncType')

% Inter requires multiple bodies
if strcmp(cfg.syncType, 'inter') && numel(varargin) <= 2,
    error('Not enough argument, Min: tsvFile1, tsvFile2, cfg')
end

% dataType - marker or feature
errorIfNotField(cfg, 'dataType')

% feature requires feature field for all bodies
if strcmp(cfg.dataType, 'feature') 
    for tsvFileIdx = 1:numel(tsvFileIdxList),
        if ~isfield(tsvFileAll{tsvFileIdx}, 'processing')
            error(['No feature calculated for tsvFile ', num2str(tsvFileIdx)]);
        end
    end
end

% dataList - list of marker or feature used in analysis
if ~isfield(cfg, 'dataList'),
    disp('cfg.dataList was not setup. Default: All feature/marker')
    if strcmp(cfg.dataType, 'feature')
        cfg.dataList = fieldnames(tsvFileAll{1}.processing);
    elseif strcmp(cfg.dataType, 'marker')
        cfg.dataList = fieldnames(tsvFileAll{1}.markerName);
    end
end

% Display progress
dispstat('','init'); %one time only init 
dispstat('Begining the processing of event sync...','keepthis','timespamp'); 
% if strcmp(cfg_main.event_sync_global, 'inter') && length(cfg_main.prefix_body_markers) == 1, 
%     error('Synchronisation across bodies is not possible if only one body')
% elseif strcmp(cfg_main.event_sync_global, 'intra')
%     if strcmp(cfg_main.event_sync_type, 'markers') && length(cfg_main.event_sync_markers) == 1
%         error('Synchronisation intra is not possible if only one marker')
%     elseif strcmp(cfg_main.event_sync_type, 'feature') && length(cfg_main.event_sync_feature) == 1
%         error('Synchronisation intra is not possible if only one feature')
%     end
% end
% 
% % Defauls
% if isempty(cfg_main.event_sync_markers)
%     cfg_main.event_sync_markers = cfg_main.markersname_noprefix;
% end
% if isempty(cfg_main.event_sync_feature)
%     cfg_main.event_sync_feature = cfg_main.feature_list;
% end




%% COMPUTATION AREA

% INTERPERSONNEL EVENT SYNC TYPE
if strcmp(cfg.syncType, 'inter')
    % Define the permutations between all the different bodies
    permutationBody = nchoosek(1:length(tsvFileIdxList),2);
    syncQ = cell(length(tsvFileIdxList),1);
    syncq = cell(length(tsvFileIdxList),1);
    for fileIdx = 1:length(tsvFileIdxList),
        % Display progress
        percentage = (fileIdx-1)/length(tsvFileIdx)*100;
        dispstat(sprintf('Processing %i%%',int64(round(percentage))));
        % Init
        if strcmp(cfg.dataType, 'marker')
            syncQ{fileIdx} = cell(length(cfg.dataList),1);
            syncq{fileIdx} = cell(length(cfg.dataList),1);
        elseif strcmp(cfg.dataType, 'feature')
            syncQ{fileIdx} = cell(length(cfg.dataList),1);
            syncq{fileIdx} = cell(length(cfg.dataList),1);
        end
        
        for pairIdx = 1:size(permutationBody,1),
            % Select two bodies
            bodyIdx1 = permutationBody(pairIdx,1);
            bodyIdx2 = permutationBody(pairIdx,2);
            % Apply Sync on either the markers or the feature
            % -------------------------------------------------
            % MARKERS
            if strcmp(cfg.dataType, 'marker')
                tsvFile1 = tsvFileAll{bodyIdx1};
                tsvFile2 = tsvFileAll{bodyIdx2};
                % For every marker assigned in cfg
                [indexList] = findIndexList(tsvFile1.markerName, cfg.dataList);
                for markerPairIdx = 1:length(cfg.dataList)
                    % Init 
                    Qvector = [];%zeros(1,length(cfg.event_sync_markers)*3);
                    qvector = [];%zeros(1,length(cfg.event_sync_markers)*3);
                    if isempty(syncQ{fileIdx}{markerPairIdx})
                        syncQ{fileIdx}{markerPairIdx} = zeros(length(tsvFileIdxList));
                        syncq{fileIdx}{markerPairIdx} = zeros(length(tsvFileIdxList));
                    end
                    % Find index in 3D
                    idx = ((indexList(markerPairIdx)-1)*3)+1;
                    % Go Through the coordinate and process sync
                    for k = idx:idx+2
                        [Qtmp, qtmp] = event_synchronization_on_peaks(tsvFile1.data(:,k), tsvFile2.data(:,k), 'tot', 1);
                        Qvector(end+1) = Qtmp;
                        qvector(end+1) = qtmp;
                    end
                    syncQ{fileIdx}{markerPairIdx}(bodyIdx1,bodyIdx2) = nanmean(Qvector);
                    syncQ{fileIdx}{markerPairIdx}(bodyIdx2,bodyIdx1) = nanmean(Qvector);
                    syncq{fileIdx}{markerPairIdx}(bodyIdx1,bodyIdx2) = nanmean(qvector);
                    syncq{fileIdx}{markerPairIdx}(bodyIdx2,bodyIdx1) = nanmean(qvector); 
                end
            % -------------------------------------------------
            % feature
            elseif strcmp(cfg.dataType, 'feature')
                % For every feature assigned in cfg
                for i = 1:length(cfg.dataList),
                    featBody1 = tsvFileAll{bodyIdx1}.processing.(cfg.dataList{i});
                    featBody2 = tsvFileAll{bodyIdx2}.processing.(cfg.dataList{i});
                    
                    % Check if it is a table from the "continuous"
                    % parameter
                    if size(featBody1,1) == 1
                        error(['Feature ' cfg.dataList{i} 'is just a single value. Use feature with continous value'])
                    elseif size(featBody1,2) > 2
                        error(['Feature ' cfg.dataList{i} 'is a multicolumns table. Use feature with vectors'])
                    end
                    % Init
                    if isempty(syncQ{fileIdx}{i})
                        syncQ{fileIdx}{i} = zeros(length(tsvFileIdxList));
                        syncq{fileIdx}{i} = zeros(length(tsvFileIdxList));
                    end
                    % Go Through the columns of the feature and process sync
                    [Qtmp, qtmp] = event_synchronization_on_peaks(featBody1, featBody2);
                    syncQ{fileIdx}{i}(bodyIdx1,bodyIdx2) = Qtmp;
                    syncQ{fileIdx}{i}(bodyIdx2,bodyIdx1) = Qtmp;
                    syncq{fileIdx}{i}(bodyIdx1,bodyIdx2) = qtmp;
                    syncq{fileIdx}{i}(bodyIdx2,bodyIdx1) = qtmp;
                end
            end
        end
        % Display Figures
        if cfg.display
            if strcmp(cfg.dataType, 'marker')
                for markerPairIdx = 1:length(cfg.dataList)
                    figure('Name',['Sync INTER Markers for file ' num2str(fileIdx) 'and marker ' cfg.dataList{markerPairIdx}],'NumberTitle','off')
                    imagesc(syncQ{fileIdx}{markerPairIdx})
                    set(gca,'XTick',1:length(tsvFileIdxList))
                    set(gca,'YTick',1:length(tsvFileIdxList))
                    set(gca,'xticklabel',tsvFileIdxList)
                    set(gca,'yticklabel',tsvFileIdxList)
                    colorbar
                end
            elseif strcmp(cfg.dataType, 'feature')
                for fatureIdx = 1:length(cfg.dataList)
                    figure('Name',['Sync INTER feature for file ' num2str(fileIdx) 'and marker ' cfg.dataList{fatureIdx}],'NumberTitle','off')
                    imagesc(syncQ{fileIdx}{fatureIdx})
                    set(gca,'XTick',1:length(tsvFileIdxList))
                    set(gca,'YTick',1:length(tsvFileIdxList))
                    set(gca,'xticklabel',tsvFileIdxList)
                    set(gca,'yticklabel',tsvFileIdxList)
                    rotateticklabel_imagesc(gca,270);
                    colorbar
                end
            end
        end
    end
    dispstat('Processing 100%%');
    %% INTRA PERSONNAL EVENT SYNC TYPE
elseif strcmp(cfg.syncType, 'intra')
    % Initialisation
    syncQ = cell(length(tsvFileIdxList),1);
    syncq = cell(length(tsvFileIdxList),1);
    for fileIdx = 1:length(tsvFileIdxList),
        % Display progress
        percentage = (fileIdx-1)/length(tsvFileIdxList)*100;
        dispstat(sprintf('Processing %i%%',int64(round(percentage))));
        
        % Init
        syncQ{fileIdx} = zeros(length(cfg.dataList));
        syncq{fileIdx} = zeros(length(cfg.dataList));
        % Define TSV file & feature
        tsvFile = tsvFileAll{fileIdx};
        % Apply Sync on either the markers or the feature
        % ------------------------------------------------------------
        % MARKERS
        if strcmp(cfg.dataType, 'marker')
            % For every marker assigned in cfg
            [indexList] = findIndexList(tsvFile.markerName, cfg.dataList);
            tsvFileMarker = mcgetmarker(tsvFile,indexList);
         
            % Define the permutations between all the different markers
            permutationMarker = nchoosek(1:length(cfg.dataList),2);
            % Computation of Event sync for intra
            for markerPairIdx = 1:size(permutationMarker,1),
                % Select two feature
                markerIdx1 = permutationMarker(markerPairIdx,1);
                markerIdx2 = permutationMarker(markerPairIdx,2);
                % Find index in 3D
                idx1 = ((markerIdx1-1)*3)+1;
                idx2 = ((markerIdx2-1)*3)+1;
                marker1 = tsvFileMarker.data(:,idx1:idx1+2);
                marker2 = tsvFileMarker.data(:,idx2:idx2+2);
                
                % Init
                Qvector = [];
                qvector = [];
                % Compute for each dimension
                for k = 1:3
                    [Qtmp, qtmp] = event_synchronization_on_peaks(marker1(:,k), marker2(:,k), 'tot', 1);
                    Qvector(end+1) = Qtmp;
                    qvector(end+1) = qtmp;
                end
                % Average
                syncQ{fileIdx}(markerIdx1,markerIdx2) = nanmean(Qvector);
                syncQ{fileIdx}(markerIdx2,markerIdx1) = nanmean(Qvector);
                syncq{fileIdx}(markerIdx1,markerIdx2) = nanmean(qvector);
                syncq{fileIdx}(markerIdx2,markerIdx1) = nanmean(qvector);
            end
            % ------------------------------------------------------------
            % feature
        elseif strcmp(cfg.dataType, 'feature')
            % Define the permutations between all the different
            % feature
            permutationFeature = nchoosek(1:length(cfg.dataList),2);
            
            % Computation of Event sync for intra
            for pairIdx = 1:size(permutationFeature,1),
                % Select two feature
                featIdx1 = permutationFeature(pairIdx,1);
                featIdx2 = permutationFeature(pairIdx,2);
                feature1 = tsvFile.processing.(cfg.dataList{featIdx1});
                feature2 = tsvFile.processing.(cfg.dataList{featIdx2});
                % Check Feature 1
                if size(feature1,1) == 1,
                    error(['Feature ' cfg.dataList{featIdx1} 'is just a single value. Use feature with continous value'])
                elseif size(feature1,2) > 1,
                    error(['Feature ' cfg.dataList{featIdx1} 'is a multicolumns table. Use feature with vectors'])
                end
                % Check Feature 2
                if size(feature2,1) == 1,
                    error(['Feature ' cfg.dataList{featIdx2} 'is just a single value. Use feature with continous value'])
                elseif size(feature2,2) > 1,
                    error(['Feature ' cfg.dataList{featIdx2} 'is a multicolumns table. Use feature with vectors'])
                end
                
                % Process the event sync
                [Qtmp, qtmp] = event_synchronization_on_peaks(feature1, feature2, 'tot', 1);
                % Place it
                syncQ{fileIdx}(featIdx1,featIdx2) = Qtmp;
                syncQ{fileIdx}(featIdx2,featIdx1) = Qtmp;
                syncq{fileIdx}(featIdx1,featIdx2) = qtmp;
                syncq{fileIdx}(featIdx2,featIdx1) = qtmp;
            end
        end
        if cfg.display
            if strcmp(cfg.dataType, 'marker')
                figure('Name',['Sync INTRA Markers for file ' tsvFile.info.filename],'NumberTitle','off')
            elseif strcmp(cfg.dataType, 'feature')
                figure('Name',['Sync INTRA Feature for file ' tsvFile.info.filename],'NumberTitle','off')               
            end
            axisValues = cfg.dataList;
            imagesc(syncQ{fileIdx})
            set(gca,'XTick',1:length(axisValues))
            set(gca,'YTick',1:length(axisValues))
            set(gca,'xticklabel',axisValues)
            set(gca,'yticklabel',axisValues)
            rotateticklabel_imagesc(gca,270);
            colorbar
        end
    end
    dispstat('Processing 100%%');
end


end

