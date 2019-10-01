%% (3) Graph Analysis and Metrics
%% Intro
% _Created by Brighton Ancelin_
%
% This toolbox has four main objectives:
% 
% # Help the user import a graph from a file or matrix into MATLAB as a graph object
% # Help the user perform common operations and alterations on graph objects
% # *Provide tools for basic analysis of graph objects and calculation of key metrics*
% # Provide a variety of different visualization means for graph objects, both time-varying and time-invariant
%
% In this page, we will discuss how to generate a graph metrics struct.
%% All-in-One Metrics Struct
% A metrics struct contains many different fields of data about a graph. A
% brief overview of each field can be seen by using the MATLAB command: 
% help getGraphMetrics
% 
% There are two main functions associated with graph metrics fields:
% *getGraphMetrics* and *exportMetricsFile*. The prior can be used to
% obtain a metrics struct, and the latter can be used to export said
% metric data into a text file. The latter function will call the prior
% internally and return the struct, so there's seldom a reason to call the
% two in succession. An example of proper usage is as follows:
metrics = exportMetricsFile(graph([1,2,3],[2,3,1]),'Title of Graph','filename');
help getGraphMetrics;