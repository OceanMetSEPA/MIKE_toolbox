function [ triIndex,nodeIndex ] = meshIndex(xp,yp,tri,xMesh,yMesh)
% Find mesh triangle(s) containing point(s), and nearest grid node
%
% INPUTS:
% xp,yp : eastings/northings of points to locate;
%   Either: tri,x,y : mesh triangles and node points OR
%           mikeMesh : struct containing above info (as returned by MIKE.loadMesh function) OR
%           dfsuStruct : matlab.io.MatFile (as returned by dfsu2Struct)
% OUTPUTS:
% triIndex: index(ices) of mesh triangles containing points
% nodeIndex: nearest node
%
% NOTE:
% any points not within mesh return NaN
%
% UPDATE (20180801)
% Rewritten to take advantage of 'tsearchn' function which is much faster
% than previous code, but returned NaNs near boundary in 2013 matlab
% version.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Workfile:   meshIndex.m  $
% $Revision:   1.2  $
% $Author:   Ted.Schlicke  $
% $Date:   Aug 01 2018 15:29:10  $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    help MIKE.meshIndex
    return
end
origShape=size(xp);
xp=xp(:);
yp=yp(:);

if nargin==3 % xmesh,ymesh,tri bundled into struct?
    fn=fieldnames(tri);
    try        
        xfield=fn{contains(fn,'xmesh','ig',1)};
        xMesh=tri.(xfield);
        yfield=fn{contains(fn,'ymesh','ig',1)};
        yMesh=tri.(yfield);
        triField=fn{contains(fn,'triMesh','ig',1)|contains(fn,'MeshIndices','ig',1)|contains(fn,'tri')};
        tri=tri.(triField);
    catch
        disp(fn)
        error('Can''t find mesh info from argument 3')
    end
end

xy=[xMesh,yMesh];
triIndex=tsearchn(xy,tri,[xp,yp]);
triIndex=reshape(triIndex,origShape);
if nargout>1
    nodeIndex=kNearestNeighbors([xMesh,yMesh],[xp,yp],1);
end

end
