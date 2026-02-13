function rootPath = setpath
%% Project Initialization Script
% Save this file in: ...\Standard_Model_Fitting\setPath.m

% 1. Get the directory where THIS script is saved (The Root)
rootPath = fileparts(mfilename('fullpath'));

% 2. Add all subfolders (Run, Data, Functions, etc.) to MATLAB's path
addpath(genpath(rootPath));

% 3. Set a global or base variable for your paths
assignin('base', 'baseDir', rootPath);

fprintf('Project root set to: %s\n', rootPath);
fprintf('All subfolders added to path.\n');