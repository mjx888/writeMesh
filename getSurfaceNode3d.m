function nodes = getSurfaceNode3d( ele )
% getSurfaceNodes: find nodes at external surface in a mesh
% Works for tetrahedral and hexahedral mesh.
% Works for linear and quadratic element.
%
% Input:
%   ele - Mesh elements. 
%
% Output:
%   nodes - column vector. Node indices for nodes at external surface in 
%           the mesh.
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%
    
    % ---------------------------------------------------------------------
    ele_wid = size( ele, 2 );

    if ele_wid == 4
        % Four faces per tet
        faces = [ ele(:,[2 3 4]);   
                  ele(:,[1 4 3]);   
                  ele(:,[1 2 4]);  
                  ele(:,[1 3 2]) ];
    elseif ele_wid == 8
        % Build all 6 faces per hex
        faces = [ ele(:,[1 2 3 4]);   % bottom
                  ele(:,[5 6 7 8]);   % top
                  ele(:,[1 2 6 5]);   % side
                  ele(:,[2 3 7 6]);   % side
                  ele(:,[3 4 8 7]);   % side
                  ele(:,[4 1 5 8]) ]; % side
    elseif ele_wid == 10
        % Four faces per tet
        faces = [ ele(:,[2 6 3 10 4 9]);   
                  ele(:,[1 8 4 10 3 7]);   
                  ele(:,[1 5 2 9 4 8]);  
                  ele(:,[1 7 3 6 2 5]) ];
    elseif ele_wid == 20
        % Build all 6 faces per hex
        faces = [ ele(:,[1 9 2 10 3 11 4 12]);   % bottom
                  ele(:,[5 13 6 14 7 15 8 16]);   % top
                  ele(:,[1 9 2 18 6 13 5 17]);   % side
                  ele(:,[2 10 3 19 7 14 6 18]);   % side
                  ele(:,[3 11 4 20 8 15 7 19]);   % side
                  ele(:,[4 12 1 17 5 16 8 20]) ]; % side
    else
        error("bad cases");
    end
    
    % ---------------------------------------------------------------------
    % Deduplicate by treating faces as unordered for matching
    key = sort(faces, 2);
    [~, ~, ic] = unique(key, 'rows', 'stable');
    counts = accumarray(ic, 1);
    isBoundary = counts(ic) == 1;
    
    % Keep boundary faces (appear exactly once)
    F = faces(isBoundary, :);
    
    nodes = unique(F(:));
    % ---------------------------------------------------------------------
end