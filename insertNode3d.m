function [vertU, eleU] = insertNode3d(vert, ele)
% insertNode3d: inserts midpoints into all edges of linear elements to form 
% quadratic elements.
%
% Works for tetrahedral and hexahedral element.
%
%   [vertU, eleU] = insertNode(vert, ele);
%
% input:
%   vert : (N x 3) array of node coordinates
%          vert(i,:) = (x, y, z) for the i-th node
%
%   ele : (M x 4) array of tet OR (M x 8) array of hex
%         ele(k,:) = [i1, i2, i3, i4] (indices of the tet's nodes)
%
% output:
%   vertU : (N + E) x 2 array of updated node coordinates
%           (original nodes + newly inserted midpoints)
%
%   eleU : (M x 10) array of quadratic tetrahedral connectivity
%          eleU(k,:) = [i1, i2, i3, i4, i5, i6, i7, i8, i9, i10],
%          where the last 6 entries are the midpoint node indices on edges
%          OR
%          (M x 20) array of quadratic hexahedral connectivity
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    ele_wid = size(ele,2);

    if ele_wid == 4
        % linear tetrahedral
        [vertU, eleU] = insertNodeTet(vert, ele);

    elseif ele_wid == 8
        % linear hexahedral
        [vertU, eleU] = insertNodeHex(vert, ele);

    else
        error('Number of columns in ele should be 4 or 8.');
    end
    
end

function [vertU, eleU] = insertNodeTet(vert, ele)
% insertNodeTet: inserts midpoints into all tetrahedral edges to form 
% 10-node quadratic elements.

    % Number of original vertices and elements
    nV = size(vert,1);
    nE = size(ele,1);
    
    % -------------------------------------------------------------
    % 1) Gather all midpoints
    % -------------------------------------------------------------
    % Vectorize by listing the edges in exactly that order:
    i1 = reshape(ele(:,[1,2,3,4,4,4])', [], 1);   % 6*nE-by-1
    i2 = reshape(ele(:,[2,3,1,1,2,3])', [], 1);   % 6*nE-by-1

    % Compute all midpoints in a single shot (still in the same order)
    newCoords = 0.5* (vert(i1, :) + vert(i2, :));  % 6*nE-by-3

    % -------------------------------------------------------------
    % 2) Combine the original vertices + newly created midpoints
    % -------------------------------------------------------------
    V2 = [vert; newCoords];  % (nV + 6*nE)-by-3

    % -------------------------------------------------------------
    % 3) Remove duplicates with "stable" to replicate original order
    %    The output 'ic' tells us, for each row in V2,
    %    which unique index it corresponds to in vertU.
    % -------------------------------------------------------------
    [vertU, ~, ic] = unique(V2, 'stable', 'rows');
    
    % The first occurrence of any coordinate is kept, in the
    % order they appear in V2.

    % -------------------------------------------------------------
    % 4) Construct the final 10-node connectivity "eleU"
    %    The first 4 columns should be the mapped indices of the
    %    original corners.  The next 6 columns are the midpoints.
    %    Then do a stable UNIQUE and update indices with "eleU==i => ic(i)"
    % -------------------------------------------------------------

    % Vector-map the original corners:
    %  If ele is nE-by-4, then ic(ele) is also nE-by-4.
    cornerMapped = ic(ele);  % each corner i => ic(i)

    % Vector-map the midpoints:
    %   The last 6*nE rows in V2 are the newly created midpoints,
    %   and their final indices in vertU are ic(nV+1 : nV+6*nE).
    %   We reshape them by 6 per element.
    midMapped = reshape(ic(nV+1:end), 6, nE).';  % nE-by-6

    % Combine corners + midpoints
    eleU = [cornerMapped, midMapped];
end

function [vertU, eleU] = insertNodeHex(vert, ele)
% insertNodeTet: inserts midpoints into all hexahedral edges to form 
% 20-node quadratic elements.

    % Number of original vertices and elements
    nV = size(vert,1);
    nE = size(ele,1);
    
    % -------------------------------------------------------------
    % 1) Gather all midpoints
    % -------------------------------------------------------------
    % Vectorize by listing the edges in exactly that order:
    i1 = reshape(ele(:,[1,2,3,4,5,6,7,8,1,2,3,4])', [], 1);   % 12*nE-by-1
    i2 = reshape(ele(:,[2,3,4,1,6,7,8,5,5,6,7,8])', [], 1);   % 12*nE-by-1
    
    % Compute all midpoints in a single shot (still in the same order)
    newCoords = 0.5* (vert(i1, :) + vert(i2, :));  % 12*nE-by-3

    % -------------------------------------------------------------
    % 2) Combine the original vertices + newly created midpoints
    % -------------------------------------------------------------
    V2 = [vert; newCoords];  % (nV + 12*nE)-by-3

    % -------------------------------------------------------------
    % 3) Remove duplicates with "stable" to replicate original order
    %    The output 'ic' tells us, for each row in V2,
    %    which unique index it corresponds to in vertU.
    % -------------------------------------------------------------
    [vertU, ~, ic] = unique(V2, 'stable', 'rows');
    
    % The first occurrence of any coordinate is kept, in the
    % order they appear in V2.

    % -------------------------------------------------------------
    % 4) Construct the final 20-node connectivity "eleU"
    %    The first 8 columns should be the mapped indices of the
    %    original corners.  The next 12 columns are the midpoints.
    %    Then do a stable UNIQUE and update indices with "eleU==i => ic(i)"
    % -------------------------------------------------------------

    % Vector-map the original corners:
    %  If ele is nE-by-8, then ic(ele) is also nE-by-8.
    cornerMapped = ic(ele);  % each corner i => ic(i)

    % Vector-map the midpoints:
    %   The last 12*nE rows in V2 are the newly created midpoints,
    %   and their final indices in vertU are ic(nV+1 : nV+12*nE).
    %   We reshape them by 12 per element.
    midMapped = reshape(ic(nV+1:end), 12, nE).';  % nE-by-12

    % Combine corners + midpoints
    eleU = [cornerMapped, midMapped];
end