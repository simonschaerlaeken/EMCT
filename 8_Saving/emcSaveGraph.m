function [] = emcSaveGraph(outputDir)
% Saves all the graphs opened and active
% 
% syntax
% emcSaveGraph(outputDir);
% 
% input parameters
% outputDir: str path where to save all graphs
% 
% output
% -
% 
% examples
% emcSaveGraph('C:/save');
% 
% comments
% Saves every graph using the name given when opening the feature.
% Saves graph in pdf format
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
% Check if exists
if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end
%% COMPUTING AREA
disp('[SAVE] Saving all available graphs into PDFs')
h = get(0,'children');
for i=1:length(h)
    figName = get(h(i),'Name');
    figPath = [outputDir, filesep, figName,'__', date, '.fig'];
    disp(['Saving..', figPath, ])
    saveas(h(i),figPath, 'fig');
end
end

