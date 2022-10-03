function [dfsuData,dfsuFile]=getDfsuInfo(dfsuFileName,varargin)
%% Get information about dfsuFile
% Input parameters:
% dfsuFileName - .dfsu file to load (HD, PART)
%
% If dfsu file generated in run using MPI, index order of x,y,z coordinates will match order of
% triangles in mesh. That is, point 'i' of dfsu ouput will be within
% triangle 'i' of mesh.
% If run used OpenMP, index order doesn't match (?)
%
% Can be useful to have these matched up, e.g. for calculating loads (conc
% * area * depth), so function can determine this order if it is passed the
% relevant mesh file
%
% The x,y,z coordinates returned by this function are in order found in
% .dfsu file.
%
% INPUT:
% .dfsu file 
%
% OUTPUT:
% struct containing coordinate / mesh information about MIKE .dfsu file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   getDfsuInfo.m  $
% $Revision:   1.0  $
% $Author:   Ted.Schlicke  $
% $Date:   Sep 25 2018 09:42:24  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%error('Use MIKE.readDfsuFile with ''data'',0. This copes with 3d dfsu files') 

if nargin==0
    help MIKE.getDfsuInfo
    return
end

% 20220310 - Changed below
% NET.addAssembly('DHI.Generic.MikeZero.DFS');
% import DHI.Generic.MikeZero.DFS.*;
NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
NETaddAssembly('DHI.Generic.MikeZero.EUM.dll');
import DHI.Generic.MikeZero.DFS.*;
import DHI.Generic.MikeZero.*.*;

options=struct;
options.meshFile=[]; % use this to determine index order
options.nodeElements=false; % return extra coordinate info from dfsu file
options=checkArguments(options,varargin);

cprintf('g','Loading ''%s''\n',dfsuFileName)
dfsuFile=DfsFileFactory.DfsuFileOpen(dfsuFileName);   % Open Mike Output file

% Spatial information:
% Node points
xn = double(dfsuFile.X);
yn = double(dfsuFile.Y);
zn = double(dfsuFile.Z);
% Create element table in Matlab format
tn = mzNetFromElmtArray(dfsuFile.ElementTable);
% also calculate element center coordinates
[xe,ye,ze] = mzCalcElmtCenterCoords(tn,xn,yn,zn);
% Triangulate grid:
tri = mzTriangulateElmtCenters(xe,ye,tn); % triangles defining grid
cprintf('cyan','%d nodes and %d triangles\n',length(xn),length(tri))
% Model parameters:
Ni=dfsuFile.ItemInfo.Count; % number of parameters in file
parameterArray=cell(Ni,1);
for i=1:Ni
    di=dfsuFile.ItemInfo.Item(i-1); % i'th item
    diq=di.Quantity;
    si=struct('Name',char(di.Name),'Index',i,'Item',char(diq.Item),'ItemDescription',char(diq.ItemDescription),'Unit',char(diq.Unit),'UnitDescription',char(diq.UnitDescription));
    parameterArray{i}=si;
end
parameterArray=vertcat(parameterArray{:});

%% Time stuff - find start time and time sep
Nt = double(dfsuFile.NumberOfTimeSteps);
%cprintf('mag','File has %d timesteps\n',Nt)
t=dfsuFile.StartDateTime;
t0=datenum(double(t.Year),double(t.Month),double(t.Day),double(t.Hour),double(t.Minute),double(t.Second)); % Start of model run
timeSep=dfsuFile.TimeStepInSeconds/(24*60*60); % Fractions of a day
modelTimes=(t0+double(0:(Nt-1))*timeSep)'; % time steps as datenums
%fprintf('Model Run Dates:\n%s to\n%s\n',datestr(min(modelTimes)),datestr(max(modelTimes)))

%% Get index order if mesh file supplied:
if ~isempty(options.meshFile)
    try
        if ischar(options.meshFile)
            fprintf('Checking mesh index order\n')
            [triMesh,nodes]=mzReadMesh(options.meshFile);
            xMesh=nodes(:,1);
            yMesh=nodes(:,2);
        elseif isstruct(options.meshFile)
            xMesh=options.meshFile.xMesh;
            yMesh=options.meshFile.yMesh;
            triMesh=options.meshFile.triMesh;
        end
    catch err
        disp(err)
        error('Problem reading mesh info :-(')
    end
    % Find index order
    indexOrder=MIKE.meshIndex(xe,ye,triMesh,xMesh,yMesh);
    di=diff(indexOrder);
    if ~all(di==1)
        warning('Index order alert!!!!!!!')
    else
        %fprintf('Index order ok\n')
    end
    [~,indexOrder]=sort(indexOrder);
else
    indexOrder=[];
end

% Now get file size
di=dir(dfsuFileName);
% Prepare struct with model info
dfsuData=struct('Name',dfsuFileName,'Time',modelTimes,'X',xe,'Y',ye,'Z',ze,'Tri',tri,'Parameters',parameterArray,'IndexOrder',indexOrder,'GigaBytes',di.bytes/1e9,'DateCreated',di.date);
if options.nodeElements
    dfsuData.XDfsu=xn(:);
    dfsuData.YDfsu=yn(:);
    dfsuData.ZDfsu=zn(:);
    dfsuData.Elements=tn;
end
dfsuFile.Close();