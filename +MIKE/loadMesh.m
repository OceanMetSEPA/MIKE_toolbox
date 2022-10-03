function [ mikeMesh ] = loadMesh( meshFileName,varargin )
% Load MIKE mesh into a struct
%
% INPUT:
% meshFileName: either a MIKE .mesh file or a directory containing a single .mesh file
%
% Optional Inputs:
% boundary: calculate mesh boundary as well (NB this can be quite slow for large meshes)
% plot: plot mesh
%
% OUTPUT:
% struct containing fields:
%   triMesh - triangulation
%   xMesh - x coordinates of nodes
%   yMesh - y coordinates of nodes
%   zMesh - z coordinates of nodes
%   codes - MIKE codes indicating whether node is internal, land, open boundary etc
%   proj  - projection of mesh
%   zunit - unit of z coordinate (1000 used to denote meters)
%   nodes - matrix of [xMesh,yMesh,zMesh,codes]
% 
% If 'boundary' option is true: 
%  boundaryIndices - indices of coordinates on boundary
%  xMeshBoundary - x coordinates of boundary nodes
%  yMeshBoundary - y coordinates of boundary nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   loadMesh.m  $
% $Revision:   1.0  $
% $Author:   Ted.Schlicke  $
% $Date:   Jul 20 2018 10:02:20  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    help MIKE.loadMesh
    return
end

options=struct;
options.plot=false;
options.boundary=true;
options.convert=false;
options=checkArguments(options,varargin);

if isdir(meshFileName)
    meshFileDir=meshFileName;
    meshFileName=fileFinder(meshFileDir,'.mesh');
    if isempty(meshFileName)
        error('No mesh files found in ''%s''',meshFileDir)
    elseif iscellstr(meshFileName)
        disp(meshFileName);
        error('Multiple mesh files found!')
    end
end

[triMesh,nodes,proj,zUnitKey]=mzReadMesh(meshFileName);
xMesh=nodes(:,1);
yMesh=nodes(:,2);
zMesh=nodes(:,3);
codes=nodes(:,4);

if options.convert
    tic
    [xMesh,yMesh]=OS.catCoordinates(xMesh,yMesh);
    toc
    nodes(:,1)=xMesh;
    nodes(:,2)=yMesh;
end

mikeMesh=struct('triMesh',triMesh,'xMesh',xMesh,'yMesh',yMesh,'zMesh',zMesh,'codes',codes,'proj',proj,'zunit',zUnitKey,'nodes',nodes);

if options.boundary
    [boundaryIndices,xMeshBoundary,yMeshBoundary]=meshBoundary(triMesh,xMesh,yMesh);
    mikeMesh.boundaryIndices=boundaryIndices;
    mikeMesh.xMeshBoundary=xMeshBoundary;
    mikeMesh.yMeshBoundary=yMeshBoundary;
end

if options.plot
    prepareFigure
    trisurf(mikeMesh.triMesh,mikeMesh.xMesh,mikeMesh.yMesh,mikeMesh.zMesh,'EdgeColor','none')
    try
        plot(xMeshBoundary,yMeshBoundary,'-k','LineWidth',2)
    catch
        % oh well
    end
end


end
