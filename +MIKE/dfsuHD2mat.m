function [ hdMatFile ] = dfsuHD2mat( dfsuFileName,matFileName,varargin)
% Generate a .mat file from a MIKE .dfsu file
% Much smaller, and easier to extract timeseries
%
% NB!!! This can be slow (several hours sometimes if dfsu size >> 1GB)
%
% Store model output in transposed matrices, as it can be much faster to
% extract columns than rows from large matrices.
%
% param          size = [Np x Nt]
% paramTimeRow : size = [Nt x Np]
%
% where Nt = number of timesteps
%       Np = number of points
%
% To extract all spatial data at single timestep, use:
% spatialDataAtSingleTimestep = param(:,timeIndex);
%
% To extract time-series of data at single point:
% pointDataTimeSeries = paramTimeRow(:,pointIndex)
%
% INPUT:
% dfsuFileName - .dfsu created from MIKE hydrodynamic run
% matFileName - .mat file where output will be stored
%
% Optional Inputs:
% ts [0] - timesteps to extract (default = 0 means extract everything)
% dt [1] - extract every dt'th timestep
% verbose (true) - give user messages about how things are going
%
% OUTPUT:
% matfile - containing coordinate and velocity data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   dfsuHD2mat.m  $
% $Revision:   1.0  $
% $Author:   Ted.Schlicke  $
% $Date:   Sep 25 2018 09:42:20  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    help MIKE.dfsuHD2mat
    return
end

% Some options:
options=struct;
options.ts=0;
options.dt=1;
options.verbose=true;
% Note - next options were useful in old version of MIKE where HD index
% order didn't match up with mesh depending on memory settings
options.meshFile=[];
options.indexOrder=[];
options=checkArguments(options,varargin);

if options.verbose
    fprintf('Getting HD from file ''%s''\n',dfsuFileName)
end
dfsuInfoStruct=MIKE.getDfsuInfo(dfsuFileName,'meshFile',options.meshFile);%
Nx=length(dfsuInfoStruct.X);

% Get index order- either sequence passed to function
if isempty(options.indexOrder)
    indexOrder=dfsuInfoStruct.IndexOrder;
else
    indexOrder=options.indexOrder;
end
if isempty(indexOrder)
    indexOrder=1:Nx;
end
dfsuInfoStruct.X=dfsuInfoStruct.X(indexOrder);
dfsuInfoStruct.Y=dfsuInfoStruct.Y(indexOrder);
dfsuInfoStruct.Z=dfsuInfoStruct.Z(indexOrder);

% Sort timings:
Nt=length(dfsuInfoStruct.Time);
if options.ts==0
    ts2Extract=1:Nt;
else
    ts2Extract=options.ts;
end
if all(diff(ts2Extract)==1)
    ts2Extract=ts2Extract(1):options.dt:ts2Extract(end);
end
if options.verbose
    fprintf('%d time steps to extract\n',length(ts2Extract))
end
dfsuInfoStruct.Time=dfsuInfoStruct.Time(ts2Extract);

%%
NET.addAssembly('DHI.Generic.MikeZero.EUM');
NET.addAssembly('DHI.Generic.MikeZero.DFS');
NETaddDfsUtil();
import DHI.Generic.MikeZero.DFS.*;
import DHI.Generic.MikeZero.DFS.dfs123.*;
import DHI.Generic.MikeZero.*;
dfsuFile=DfsFileFactory.DfsuFileOpen(dfsuFileName);

fn=fieldnames(dfsuInfoStruct);
fn=stringFinder(fn,'*','nand',{'TimeIndices','Param'});

% Generate time-invariant parameters:
hdMatFile=matfile(matFileName,'Writable',true);
for i=1:length(fn)
    fni=fn{i};
    if options.verbose
        cprintf('mag','Assigning field %s\n',fni);
    end
    hdMatFile.(fni)=dfsuInfoStruct.(fni);
end

% allocate fields for time-varying parameters
%% Identify parameters we want to keep
% (everything except for 'Mass' fields- we'll store these as sparse
% matrices since there are so many zero values)
parameterStruct=dfsuInfoStruct.Parameters;
%Np=length(parameterStruct);
k=stringFinder({parameterStruct.Name},'Mass','output','bool');
parameterStruct(k)=[];
Np=length(parameterStruct);

%% Create fieldname from parameter name
parameterNames={parameterStruct.Name};
structParameterNames=cellfun(@MIKE.getFieldName,parameterNames,'Unif',0)';
% Add these fields to HD struct
for parameterIndex=1:Np
    hdMatFile.(structParameterNames{parameterIndex})=double([]);
end
hdMatFile.UVelocityTimeRow=double([]);
hdMatFile.VVelocityTimeRow=double([]);
hdMatFile.SurfaceElevationTimeRow=double([]);

%%
if options.verbose
    fprintf('Preparing indices for parameterNames:\n')
    disp(parameterNames)
end
uIndex=parameterStruct(stringFinder(parameterNames,{'U','velocity'},'output','index')).Index;
vIndex=parameterStruct(stringFinder(parameterNames,{'V','velocity'},'output','index')).Index;
seIndex=parameterStruct(stringFinder(parameterNames,'Surface','output','index')).Index;

%% Loop through timesteps getting parameter data, and updating struct
Nt=length(dfsuInfoStruct.Time);
Nx=length(dfsuInfoStruct.X);

for timeIndex=1:length(ts2Extract)
    if mod(timeIndex,100)==0 && options.verbose
        fprintf('Timestep %d of %d\n',timeIndex,Nt);
    end
    dfsuTimeStep=ts2Extract(timeIndex);
    % Get Surface Elevation:
    timeIndexData=double(dfsuFile.ReadItemTimeStep(seIndex,dfsuTimeStep-1).Data);
    timeIndexData=timeIndexData(indexOrder);
    %    hdMatFile.SurfaceElevationTimeRow(timeIndex,1:Nx)=timeIndexData;
    hdMatFile.SurfaceElevation(1:Nx,timeIndex)=timeIndexData';
    % Get UVelocity:
    timeIndexData=double(dfsuFile.ReadItemTimeStep(uIndex,dfsuTimeStep-1).Data);
    timeIndexData=timeIndexData(indexOrder);
    %    hdMatFile.UVelocityTimeRow(timeIndex,1:Nx)=timeIndexData;
    hdMatFile.UVelocity(1:Nx,timeIndex)=timeIndexData';
    % Get VVelocity:
    timeIndexData=double(dfsuFile.ReadItemTimeStep(vIndex,dfsuTimeStep-1).Data);
    timeIndexData=timeIndexData(indexOrder);
    %    hdMatFile.VVelocityTimeRow(timeIndex,1:Nx)=timeIndexData;
    hdMatFile.VVelocity(1:Nx,timeIndex)=timeIndexData';
end
if options.verbose
    fprintf('Creating transposed fields...\n')
end
hdMatFile.SurfaceElevationTimeRow=transpose(hdMatFile.SurfaceElevation);
hdMatFile.VVelocityTimeRow=transpose(hdMatFile.VVelocity);
hdMatFile.UVelocityTimeRow=transpose(hdMatFile.UVelocity);
if options.verbose
    fprintf('DONE!\n')
end

return
