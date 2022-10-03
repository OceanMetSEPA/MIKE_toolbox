% Script to load various MIKE tools for writing .dfs0 etc files
% These need to be in a script rather than in a function for some reason
%
% INPUTS: None
% OUTPUTS: None
%
% Usage: call this script from function attempting to write .dfs0, .dfsu
% etc files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   mikeImportStuff.m  $
% $Revision:   1.0  $
% $Author:   ted.schlicke  $
% $Date:   Mar 09 2018 11:26:04  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NET.addAssembly('DHI.Generic.MikeZero.EUM');
NET.addAssembly('DHI.Generic.MikeZero.DFS');
NETaddDfsUtil();
import DHI.Generic.MikeZero.DFS.*;
import DHI.Generic.MikeZero.DFS.dfs123.*;
import DHI.Generic.MikeZero.DFS.dfsu.*;
import DHI.Generic.MikeZero.*;