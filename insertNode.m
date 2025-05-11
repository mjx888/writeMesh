function [vertU, triaU] = insertNode(vert, tria)
% insertNode: inserts midpoints into all edges to form quadratic elements.
% Works for triangular and quadrilateral element
%
%   [vertU, triaU] = insertNode(vert, tria);
%
% input:
%   vert : (N x 2) array of node coordinates
%          vert(i,:) = (x, y) for the i-th node
%
%   tria : (M x 3) array of triangles OR (M x 4) array of quadrilaterals
%          tria(k,:) = [i1, i2, i3] (1-based indices of the triangle's nodes)
%
% output:
%   vertU : (N + E) x 2 array of updated node coordinates
%           (original nodes + newly inserted midpoints)
%
%   triaU : (M x 6) array of "quadratic" triangle connectivity
%           triaU(k,:) = [i1, i2, i3, i4, i5, i6],
%           where the last 3 entries are the midpoint node indices on edges
%           OR
%           (M x 8) array of "quadratic" quadrilateral connectivity
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    ele_wid = size(tria,2);

    if ele_wid == 3
        % triangle
        [vertU, triaU] = insertNodeTria(vert, tria);
    elseif ele_wid == 4
        % quadrilateral
        [vertU, triaU] = insertNodeQuad(vert, tria);
    else
        error('Wierd case of ele_wid.');
    end
    
end

function [vertU, triaU] = insertNodeTria(vert, tria)
% insertNodeTria: inserts midpoints into all triangle edges to form 6-node elements.

    % Number of original vertices and elements
    nV = size(vert,1);
    nT = size(tria,1);

    % -------------------------------------------------------------
    % 1) "tria_app" for the edges (1->2, 2->3, 3->1), same as original
    % -------------------------------------------------------------
    tria_app = [tria, tria(:,1)];  % M-by-4

    % -------------------------------------------------------------
    % 2) Gather *all* midpoints in the same j,k order
    %    The original loop was:
    %      for j=1:nT
    %         for k=1:3
    %             i1 = tria_app(j,k)
    %             i2 = tria_app(j,k+1)
    %             newVert = 0.5*(vert(i1,:) + vert(i2,:))
    % -------------------------------------------------------------
    % Vectorize this by listing the edges in exactly that order:
    i1 = reshape(tria_app(:,1:3)', [], 1);   % 3*nT-by-1
    i2 = reshape(tria_app(:,2:4)', [], 1);   % 3*nT-by-1

    % Compute all midpoints in a single shot (still in the same order)
    newCoords = 0.5* (vert(i1, :) + vert(i2, :));  % 3*nT-by-2

    % -------------------------------------------------------------
    % 3) Combine the original vertices + newly created midpoints
    % -------------------------------------------------------------
    % In your original code, you appended each midpoint to 'vert'
    % and then called UNIQUE at the very end. 
    %
    % Here we do it in one block:
    V2 = [vert; newCoords];  % (nV + 3*nT)-by-2

    % -------------------------------------------------------------
    % 4) Remove duplicates with "stable" to replicate original order
    %    The output 'ic' tells us, for each row in V2,
    %    which unique index it corresponds to in vertU.
    % -------------------------------------------------------------
    [vertU, ~, ic] = unique(V2, 'stable', 'rows');
    %
    % Now 'vertU' is exactly the final vertex list your original code would
    % have produced. The first occurrence of any coordinate is kept, in the
    % order they appear in V2.

    % -------------------------------------------------------------
    % 5) Construct the final 6-node connectivity "triaU"
    %
    %    The first 3 columns should be the mapped indices of the
    %    original corners.  The next 3 columns are the midpoints.
    %
    %    In your old code, you do something like:
    %      tria(j,4) = index of midpoint on edge(1,2)
    %      tria(j,5) = index of midpoint on edge(2,3)
    %      tria(j,6) = index of midpoint on edge(3,1)
    %
    %    Then you do a stable UNIQUE and update indices with "triaU==i => ic(i)"
    % -------------------------------------------------------------

    % Vector-map the original corners:
    %  (If tria is M-by-3, then ic(tria) is also M-by-3.)
    cornerMapped = ic(tria);  % each corner i => ic(i)

    % Vector-map the midpoints:
    %   The last 3*nT rows in V2 are the newly created midpoints,
    %   and their final indices in vertU are ic(nV+1 : nV+3*nT).
    %   We reshape them by triangle, 3 per triangle.
    midMapped = reshape(ic(nV+1:end), 3, nT).';  % M-by-3

    % Combine corners + midpoints
    triaU = [cornerMapped, midMapped];
end


function [vertU, triaU] = insertNodeQuad(vert, tria)
% insertNodeQuad: inserts midpoints into all quadrilateral edges to form 8-node elements.

    % Number of original vertices and elements
    nV = size(vert,1);
    nT = size(tria,1);

    % 1) "tria_app" for the edges (1->2, 2->3, 3->1), same as original
    tria_app = [tria, tria(:,1)];

    % Vectorize this by listing the edges in exactly that order:
    i1 = reshape(tria_app(:,1:4)', [], 1); 
    i2 = reshape(tria_app(:,2:5)', [], 1);  

    % Compute all midpoints in a single shot (still in the same order)
    newCoords = 0.5* (vert(i1, :) + vert(i2, :)); 

    % 3) Combine the original vertices + newly created midpoints
    V2 = [vert; newCoords]; 

    % 4) Remove duplicates with "stable" to replicate original order
    [vertU, ~, ic] = unique(V2, 'stable', 'rows');

    % 5) Construct the final 6-node connectivity "triaU"
    % Vector-map the original corners:
    cornerMapped = ic(tria);  % each corner i => ic(i)

    % Vector-map the midpoints:
    midMapped = reshape(ic(nV+1:end), 4, nT).';

    % Combine corners + midpoints
    triaU = [cornerMapped, midMapped];
    
end