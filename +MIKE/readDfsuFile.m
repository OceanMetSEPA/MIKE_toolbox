function [dfsuData]=readDfsuFile(dfsuFileName,varargin)
% Get spatial output from a MIKE .dfsu file
%
% Spatial output from MIKE model runs is stored in .dfsu files. This includes
% coordinate information and model predictions for the items specified in
% the .m21fm / .m3fm. This function attempts to extract this data for
% plotting / further analysis.
%
% WARNING: These files can be very large - 10s of Gigabytes or more -
% and you may run out of  memory and/or freeze your computer if you attempt
% to load the entire dataset at once. It may be worth trying to extract
% subsections of the output to start with using the optional inputs.
%
% INPUT:
%   dfsuFileName - .dfsu file to load (HD, PART)
%
% Optional Inputs:
% 'parameter' [] - parameters to return. If empty, function will return all parameters
% 'layer' [] - layers (of 3d model) to return. 1 corresponds to bottom layer
% 'index' [] - spatial indices to return. If empty, return entire domain
% 't0' [0] - return time steps after this value (either a datenum or index)
% 't1' [inf] - return timesteps before this value
% 'dt' - time separation in minutes (files can be huge so you might not want every step)
% 'squeeze' [true] - apply 'squeeze' function to data
% 'data' [true] - return data; otherwise return info about model output
%
% OUTPUT:
% 'dfsuData' - struct with fields:
%   'Name' - dfsu file name
%   'DateTime' - timesteps extracted from model
%   'X' - X coordinate of data point (NOutputElements,1)
%   'Y' - Y coordinate of data point (NOutputElements,1)
%   'Z' - Z coordinate of data point (NoutputElements,NLayers)
%   'ElementTable' - table of indices relating nodes to elements (see below)
%   'Tri' - triangulation of data points
%   'MeshIndices' - indices of model mesh
%   'XMesh' - X coordinate of model mesh nodes
%   'YMesh' - Y coordinate of model mesh nodes
%   'ZMesh' - Z coordinate of model mesh nodes
%   'ParameterInfo' - struct (array) with parameter names & units
%   Data - matrices of data with names based on item name
% _________________________________________________________________________
% Notes:
%
% If the mesh used in the model run has NMeshNodes nodes and NMeshElements
% elements (triangles/quadrilaterals), the output file will contain:
% i) NOutputNodes = (NLayers+1)*NMeshNodes % nodes
% ii)NOutputElements =  NLayers * NMeshElements % elements
% where NLayers is the number of layers in the model (1 for 2d)
%
% For 3d models, the X,Y coordinates at each element are the same for each
% layer.
%
% Most output items (Surface Elevation, Current Speed) etc have a single
% value per output element per timestep i.e. NTimeSteps * NOutputElements.
% Some items (e.g. 'Z Coordinate') have one value per output node.
%
% Output data for items with a single value per element is stored in a 3d
% matrix whose dimensions correspond to 1) spatial index 2) layer 3) timestep.
% The optional argument 'squeeze' (default = true) applies matlab's squeeze
% function to remove singleton dimensions.
%
% The dfsui file contains an 'ElementTable' which defines how the node
% indices are connected into elements. This table has NOutputElement rows;
% the number of columns depends on whether the model is 2d or 3d, and
% whether there are quadrilateral elements present or just triangles:
% i) 2d, triangles only : 3 columns
% ii) 2d, quadrilaterals present: 4 columns
% iii) 3d, triangles only : 6 columns
% iv) 3d, quadrilaterals present: 8 columns
%
% For i) above, the ElementTable is simply the triangulation as returned by
% mzReadMesh. For iii) above, the ElementTable defines how the various
% layers are connected.
%
% This function calculates a 2d ElementTable corresponding to the model
% mesh to help with plotting. The mesh nodes are also calculated.
% For plotting the output, it may be better to use the model mesh rather than
% the output coordinates. For grids consisting of triangles alone, Matlab's
% trimesh/trisurf function can be used; for mixtures of triangles &
% quadrilaterals, the patch function can be used. See examples below.
%
% %% Sample plot : triangles only
% dfsuData=solwayData{2};
% x=dfsuData.XMesh;
% y=dfsuData.YMesh;
% z=dfsuData.ZMesh;
% tri=dfsuData.MeshIndices;
% timeStep2Plot=1;
% param2Plot=dfsuData.ConcentrationComponent10(:,timeStep2Plot);
% trisurf(tri,x,y,z,'FaceVertexCData',param2Plot,'EdgeColor','none')
%
% %% Sample plot : mixture of triangles and quadrilaterals:
% meshIndices=dfsuData.MeshIndices'; % note the transpose
% NVertices=size(meshIndices,1); % 3 or 4
% zeroIndices=find(meshIndices==0); % 0 indices indicate triangles
% meshIndices(zeroIndices)=meshIndices(zeroIndices-1); % duplicate previous point
% x=dfsuData.XMesh(meshIndices);
% y=dfsuData.YMesh(meshIndices);
% z=dfsuData.ZMesh(meshIndices);
% timeStep2Plot=1;
% param2Plot=dfsuData.ConcentrationComponent10(:,timeStep2Plot);
% % reshape for patch function call
% param2Plot=repmat(param2Plot',NVertices,1);
% patch(x,y,param2Plot,'EdgeColor','none')
%
% DHI function mzPlot can also plot model output.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   readDfsuFile.m  $
% $Revision:   1.2  $
% $Author:   Ted.Schlicke  $
% $Date:   Apr 05 2019 11:10:50  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    help MIKE.readDfsuFile
    return
end

%% Process inputs
options=struct;
options.dt=[];
options.ts=0;
options.t0=0;
options.t1=Inf;
options.parameter=[];
options.Nmax=[];
options.index=[];
options.layer=[];
options.data=true;
options.squeeze=true;
options.verbose=false;
options.reshapeFields=true;

% Read optional inputs:
options=checkArguments(options,varargin);
if options.verbose
    fprintf('Loading ''%s''\n',dfsuFileName)
end

%% Need these lines for MIKE tools to work:
% OLD CODE:
%NET.addAssembly('DHI.Generic.MikeZero.DFS');
%import DHI.Generic.MikeZero.DFS.*;
% NEW CODE (10/05/2021):
NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
NETaddAssembly('DHI.Generic.MikeZero.DFS.dll');
NETaddAssembly('DHI.Generic.MikeZero.EUM.dll');
import DHI.Generic.MikeZero.DFS.*;
import DHI.Generic.MikeZero.*.*;

% Load dfsu file and establish key properties
dfsu=DfsFileFactory.DfsuFileOpen(dfsuFileName);   % Open Mike Output file

% For input mesh with NMeshNodes (points) and NMeshElements (triangles):
NOutputElements=dfsu.NumberOfElements; % NLayers*NMeshElements
NOutputNodes=dfsu.NumberOfNodes; % (NLayers+1)*NMeshNodes
NLayers = dfsu.NumberOfLayers; % 0 for 2d dfsu outpuf file
if NLayers==0
    NLayers=1;
    nodeIndices=1:NOutputNodes;
    elementIndices=1:NOutputElements;
else
    nodeIndices=1:NLayers+1:NOutputNodes;
    elementIndices=1:NLayers+1:NOutputElements;
end
NMeshElements=NOutputElements/NLayers;
NMeshNodes=NOutputNodes/(NLayers+1);

% Spatial information:
% Node points: these are the vertices of the model mesh.
% number of nodes = NMeshNodes * (NLayers + 1)
xn = double(dfsu.X)';
yn = double(dfsu.Y)';
zn = double(dfsu.Z)';
% fprintf('%d nodes found\n',NOutputNodes)
% fprintf('%d layers\n',NLayers);
% fprintf('%d elements\n',NOutputElements)
% fprintf('(%d mesh elements)\n',NMeshElements)

% Create element table in Matlab format
elementTable = mzNetFromElmtArray(dfsu.ElementTable);
% also calculate element center coordinates - these correspond to the
% position of the model output values (centres of cells)
[xElement,yElement,zElement] = mzCalcElmtCenterCoords(elementTable,xn',yn',zn');

% Get 2d mesh of output elements
% We want to extract:
% 1) mesh used in model (useful for plotting)
% 2) mesh defining connection of output nodes
%fprintf('Determining 2d mesh...\n')
NColumns=size(elementTable,2);
tri2d=elementTable;
if NColumns>4 % 3d
    tri2d(:,NColumns/2+1:end)=[];
end
k=any(ismember(tri2d,elementIndices),2);
tri2d=tri2d(k,:);
% Remap indices in 3d 'triangulation'. Only need to do this if number of
% layers >1
% ('changem' function can do this but needs mapping toolbox (?!). But it's not tricky)
if NLayers>1
    tri2dOrig=tri2d; % Copy of triangulation
    for timeIndex=1:length(elementIndices)
        tri2d(tri2dOrig==elementIndices(timeIndex))=timeIndex;
    end
end
% NB: tri2d is equal to triMesh (i.e. model mesh triangulation)
% Now determine mesh nodes:
xn2d=xn(nodeIndices);
yn2d=yn(nodeIndices);
zn2d=zn(nodeIndices);

% Reshape coordinates based on layers
xe=reshape(xElement,NLayers,[])';
ye=reshape(yElement,NLayers,[])';
ze=reshape(zElement,NLayers,[])';

%% Process Parameters

NParameters=dfsu.ItemInfo.Count; % number of parameters in file
if NParameters==0
    warning('No items found in dfsu file! Returning...')
    dfsuData=[];
    return;
end
parameterInfo=cell(NParameters,1);
for parameterIndex=1:NParameters
    parameterItem=dfsu.ItemInfo.Item(parameterIndex-1);
    parameterName=char(parameterItem.Name);
    if isempty(parameterName)
        parameterName=sprintf('Item%d',parameterIndex);
    end
    quantity=parameterItem.Quantity;
    parameterStruct=struct('Name',parameterName,'Index',parameterIndex,'Item',char(quantity.Item),'ItemDescription',char(quantity.ItemDescription),'Unit',char(quantity.Unit),'UnitDescription',char(quantity.UnitDescription));
    parameterInfo{parameterIndex}=parameterStruct;
end
parameterInfo=vertcat(parameterInfo{:});
dfsuParameterNames={parameterInfo.Name};

if ~isempty(options.parameter)
    if options.verbose
        fprintf('Filtering parameter names.....\nParameters found:\n')
        disp(dfsuParameterNames)
    end
    %    dfsuParameterNames2Extract=stringFinder(dfsuParameterNames,options.parameter,'type','or');
    dfsuParameterNames2Extract=closestStringMatch(dfsuParameterNames,options.parameter);
else
    dfsuParameterNames2Extract=dfsuParameterNames;
end
if isempty(dfsuParameterNames2Extract)
    fprintf('!!!!!\n')
    fprintf('Valid parameters are:\n')
    disp(dfsuParameterNames);
    fprintf('You requested:\n')
    disp(options.parameter)
    error('Didn''t find any of the above parameters')
end
if options.verbose
    fprintf('Loading parameters:\n')
    disp(dfsuParameterNames2Extract)
end
itemFieldNames=MIKE.getFieldName(dfsuParameterNames2Extract);
if ischar(itemFieldNames)
    itemFieldNames=cellstr(itemFieldNames);
end

%% layers to extract
if ~isempty(options.layer)
    layers2Keep=ismember(1:NLayers,options.layer);
else
    layers2Keep=true(NLayers,1);
end
if sum(layers2Keep)==0
    error('No valid layers specified; returning')
end

%% Model items to extract:
if ~isnumeric(options.t0)
    options.t0=datenumVariableFormat(options.t0);
end
if ~isnumeric(options.t1)
    options.t1=datenumVariableFormat(options.t1);
end
% Time stuff - find start time and time sep
Nt = double(dfsu.NumberOfTimeSteps);
if options.verbose
    fprintf('File has %d timesteps\n',Nt)
end
t=dfsu.StartDateTime;
t0=datenum(double(t.Year),double(t.Month),double(t.Day),double(t.Hour),double(t.Minute),double(t.Second)); % Start of model run
timeSep=dfsu.TimeStepInSeconds/(24*60*60); % Fractions of a day
modelDatenums=(t0+double(0:(Nt-1))*timeSep); % time steps
if options.verbose
    fprintf('Model Run Dates:\n%s to\n%s\n',datestr(min(modelDatenums)),datestr(max(modelDatenums)))
end
% Might want to filter output times
timeStep=unique(round(diff(modelDatenums)*24*60*60)/60); % Time step in minutes
if length(timeStep)>1
    error('Please check timesteps (uneven separation?)')
end
if options.verbose
    fprintf('Model output file timestep = %.1f minutes\n',timeStep)
end
% Define indices for this time sequence
if options.ts==0 % Default - extract everything (0 used for this purpose in D3D functions)
    if isempty(options.dt)
        dk=1; % default - extract every timestep
    else
        dk=round(options.dt/timeStep);
        if dk<1
            error('Can''t extract that time resolution! Minimum timestep = %.1f',timeStep)
        end
    end
    if options.verbose
        fprintf('Keeping every %d timestep (total = %d)...\n',dk,Nt)
        fprintf('Number of modelTimes = %d\n',length(modelDatenums))
    end
    timeIndices2Extract=find(modelDatenums>=options.t0 & modelDatenums<options.t1 & mod(0:(Nt-1),dk)==0);
else
    timeIndices2Extract=options.ts;
    if timeIndices2Extract==-1
        fprintf('Extracting final timestep %d\n',Nt)
        timeIndices2Extract=Nt;
    end
end
% Filter start time too:
k=timeIndices2Extract>length(modelDatenums);
if any(k)
    warning('Index(ices) requested exceed length of model run (%d); cropping',length(modelDatenums))
    timeIndices2Extract(k)=[];
end
if isempty(timeIndices2Extract)
    error('No valid indices found')
end
% Filter times based on periods requested:
modelDatenums=modelDatenums(timeIndices2Extract);

% Has user specified maximum number of timesteps?
if ~isempty(options.Nmax)
    if length(modelDatenums)>options.Nmax
        timeIndices2Extract=timeIndices2Extract(1:options.Nmax);
        modelDatenums=modelDatenums(1:options.Nmax);
    else
        options.Nmax=length(modelDatenums);
        if options.verbose
            fprintf('Keeping %d timesteps...\n',length(modelDatenums))
        end
    end
end
if options.verbose
    fprintf('%d model times in file\n',length(modelDatenums))
    fprintf('Extracting %d timesteps out of %d\n',length(timeIndices2Extract),length(modelDatenums))
end

%% Has user requested specific spatial indices?
index=options.index;
extractAllIndices=isempty(index);
if extractAllIndices
    indices2Extract=1:NMeshElements;
else
    indices2Extract=index;
end
x=xe(indices2Extract,:);
y=ye(indices2Extract,:);
z=ze(indices2Extract,:);

%% Prepare coordinates of output:
% x,y don't change with depth, so only keep 1st layer:
x=x(:,1);
y=y(:,1);
% filter layers:
z=z(:,layers2Keep);

%% Prepare output struct
dfsuData=struct('Name',dfsuFileName,'DateTime',modelDatenums(:),'X',x(:),'Y',y(:),'Z',z,'ElementTable',elementTable,'MeshIndices',tri2d,'XMesh',xn2d,'YMesh',yn2d,'ZMesh',zn2d,'ItemInfo',parameterInfo);

%% Load data
if options.data
    Np=length(dfsuParameterNames2Extract);
    for parameterIndex=1:Np
        p=dfsuParameterNames2Extract{parameterIndex};
        % 20170922 - in previous version, simply used parameterIndex in
        % ReadItemTimeStep function call below- BUG! Causes problems when
        % we're not extracting all parameters: instead we need to find the
        % index of the parameter we're after and use that.
        pin=stringFinder(dfsuParameterNames,p,'output','index','type','exact');
        if options.verbose
            fprintf('Loading parameter %s (id = %d)\n',p,pin)
        end
        parameterData=arrayfun(@(i)double(dfsu.ReadItemTimeStep(pin,i-1).Data)',timeIndices2Extract,'Unif',0);
        parameterData=horzcat(parameterData{:}); % 2d array, Nt columns, one row per grid point (centre of each mesh triangle)
        cprintf('green','N Output values = %d\n',NOutputElements)
        if options.reshapeFields
            if size(parameterData,1)==NOutputElements
                %            fprintf('Reshaping and permuting... %s\n',p)
                parameterData=reshape(parameterData,NLayers,NMeshElements,[]);
            elseif size(parameterData,1)==NOutputNodes
                parameterData=reshape(parameterData,NLayers+1,NMeshNodes,[]);
            end
            parameterData=permute(parameterData,[2,1,3]);
        end
        % Filter layers:
        if size(parameterData,2)==NLayers
            parameterData=parameterData(:,layers2Keep,:);
        end
        % Filter indices:
        if ~extractAllIndices && size(parameterData,1)==NMeshElements
            parameterData=parameterData(indices2Extract,:,:);
        end
        fn=itemFieldNames{parameterIndex};
        dfsuData.(fn)=parameterData;
    end
    % Squeeze?
    if options.squeeze
        Nf=length(itemFieldNames);
        for fieldIndex=1:Nf
            fni=itemFieldNames{fieldIndex};
            si=dfsuData.(fni);
            si=squeeze(si);
            dfsuData.(fni)=si;
        end
    end
end
dfsu.Close();

% All done!