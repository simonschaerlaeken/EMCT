%% (1) Importing Graphs from Files
%% Intro
% _Created by Brighton Ancelin_
% 
% This toolbox has four main objectives:
% 
% # *Help the user import a graph from a file or matrix into MATLAB as a graph object*
% # Help the user perform common operations and alterations on graph objects
% # Provide tools for basic analysis of graph objects and calculation of key metrics
% # Provide a variety of different visualization means for graph objects, both time-varying and time-invariant
%
% To get started with this toolbox, we will first focus on importing graphs 
% from files and matrices. 
%
% _Note: A filename ending in .net is assumed to be a Pajek file. Please be
% careful to change extensions if needed. More info on the Pajek format
% here:_
% <http://vlado.fmf.uni-lj.si/pub/networks/pajek/>
%% Simple Edge Files
% Graph data comes in all kinds of formats, but most commonly is the
% edge-by-row format. In this format, every line of the file contains two
% positive integers. Each line in the file represents an edge in the graph,
% and the two integers on said line are the nodes which that edge connects.
% The first integer is the "from node" and the second integer is the "to
% node." In undirected graphs, the order of these nodes does not matter, as
% edges do not have direction. In directed graphs, however, the order does
% matter, and the edge points from the "from node" to the "to node." An
% example of one such file might be the following 'ex_edges-by-rows.txt':
%
% <html>
% 1 5<br>
% 2 3<br>
% 4 6<br>
% 12 8<br>
% 9 12<br>
% </html>
% 
% Such files can easily be imported as graph objects using the *importNet*
% function. For more information on how this function works, use the MATLAB
% command: help importNet

graphObj = importNet('ex_edges-by-rows.txt',false);
graphObj.Edges
%% More Complicated Edge Files
% While many edge files have each row representing a single edge, some
% files will have more than two columns. Two of the columns will contain
% node IDs, but other data will be written as well. Some of these columns
% may be useful to import into our graph objects, such as edge weights for
% weighted graphs. For example, let's look at 'ex_complicatedEdgeFile.txt',
% whose first column contains edge weights and whose second and third
% columns contain node IDs:
%
% <html>
% 4 1 4<br>
% 2 2 3<br>
% 4 3 1<br>
% 7 5 8<br>
% 3 1 5<br>
% </html>

graphObj = importNet('ex_complicatedEdgeFile.txt',false,'nodeCols',[2,3],'weightCol',1);
graphObj.Edges
%% Duplicate Edge Handling
% Graph objects in MATLAB can't hold multiple identical edges, so
% *importNet* has four predefined methods for dealing with duplicate edges:
% 'first', 'ignore', 'sum', and 'average'. For more information on what 
% each does, use the MATLAB command: help importNet
%
% For directed graphs, the order of nodes matters, but for undirected 
% graphs it does not. Let's observe how we can import the following graph 
% using averaged weights from the file 'ex_duplicateEdges.txt':
%
% <html>
% 1 2 500<br>
% 3 4 400<br>
% 2 4 600<br>
% 2 1 200<br>
% 3 4 700<br>
% </html>

graphObj = importNet('ex_duplicateEdges.txt',false,'nodeCols',[1,2],...
		'weightCol',3,'dupMode','average');
graphObj.Edges
%% Preprocessing a Matrix for use with importNet
% If the capabilities of *importNet* just aren't cutting it, you can always
% do some preprocessing and pass *importNet* the edge matrix. For example, 
% let's say you have a file with 3 columns: the first two are node IDs, and
% the third represents the classification of the edge. We want to create a 
% graph with only classification 7 edges. Here's the file 
% 'ex_edgesClassified.txt':
%
% <html>
% 1 2 7<br>
% 3 2 3<br>
% 3 2 7<br>
% 5 7 4<br>
% 4 2 7<br>
% 2 5 7<br>
% 2 2 2<br>
% </html>
%
% Here's the code:
fileMatrix = str2num(fileread('ex_edgesClassified.txt'));
class7mask = 7 == fileMatrix(:,3);
graphObj = importNet(fileMatrix(class7mask,1:2),false);
graphObj.Edges
%% Importing a Graph as Frames
% Graphs that change over time, or temporal graphs, can be imported as a
% cell vector of graph objects. However, this is not accomplished with the
% *importNet* function, rather it is accomplished with the *getGraphFrames*
% function. This function will create graph frames from a matrix of edges
% and a column vector of timestamps corresponding to the entries in the
% edge matrix. For more information on how to use this function, use the
% MATLAB command: help getGraphFrames
data = [1,2,  0;...
	    1,2,  1;...
	    2,3,  1;...
	    4,2,  2;...
	    4,2,  3;...
	    3,2,  0;...
	    1,4,  5];
edgeMat = data(:,1:2);
timestampVec = data(:,3);
graphFrames = getGraphFrames(edgeMat,timestampVec,true,1,99);
fprintf('Frame Count: %d\n',length(graphFrames));
fprintf('Frame 1:\n'); disp(graphFrames{1}.Edges.EndNodes);
fprintf('Frame 2:\n'); disp(graphFrames{2}.Edges.EndNodes);
fprintf('Frame 3:\n'); disp(graphFrames{3}.Edges.EndNodes);
fprintf('Frame 4:\n'); disp(graphFrames{4}.Edges.EndNodes);
fprintf('Frame 5:\n'); disp(graphFrames{5}.Edges.EndNodes);
fprintf('Frame 6:\n'); disp(graphFrames{6}.Edges.EndNodes);