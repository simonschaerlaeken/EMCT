function [ tsvFile ] = emcNorm( tsvFile, cfg )
% Computes the norm and stores it as a feature 
% 
% syntax
% tsvFile = emcNorm( tsvFile, cfg )
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.-
%     [OPTIONAL]
%     *.-
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% tsvFile = emcNorm(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.norm
% 
% see also
% -
% 
%% CHECKING
% featMarker - markers on which the feature is calculated
if ~isfield(cfg, 'featMarker')||isempty(cfg.featMarker) 
    disp('cfg.featMarker was not setup. Default: all markers')
    cfg.featMarker = tsvFile.markerName;
end
%% PROCESSING
tsvTmp = mcnorm(tsvFile);
tsvFile.processing.norm = tsvTmp.data;
% PLOT
if cfg.display
    colorset = distinguishable_colors(length(cfg.featMarker));
    figure('Name',[tsvFile.info.filename, ' - Norm'],'NumberTitle','off')
    hold on
    for i = 1:size(tsvTmp.data,2)
        if strcmp(cfg.displayUnit, 's')
            time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
            plot(time, tsvTmp.data(:,i),'color', colorset(i,:))
        else
            plot(tsvTmp.data(:,i),'color', colorset(i,:))
        end
    end
    hold off
    legend(cfg.featMarker)
end
end

