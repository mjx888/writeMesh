function [ nodecoor_list, nodecoor_cell, ele_cell ] = getNodeEle3d( vert, ele, tnum )
% getNodeEle3d: get node coordinates and elements from 3d mesh
%
% Detail: Add node numbering and element numbering, and organize elements 
% into cell array. eleC{i} represent elements in the i-th phase.
%
% Works for tetrahedral and hexahedral mesh.
% Works for linear and quadratic element.
%
% usage:
%   [ nodecoor_list, nodecoor_cell, ele_cell ] = getNodeEle3d( vert, ele, tnum );
%
% input:
%   vert - Node data. N-by-3 array. x, y, z coordinates.
%       vert(i,1:3) = [x, y, z] of the i-th node
%
%   ele - Node numbering for each element
%
%   tnum - Label of material phase. P-by-1 array.
%       tnum(j,1) = k; means the j-th element is belong to the k-th phase
%
% output:
%   nodecoor_list is node coordinates, N-by-4 array.
%       nodecoor_list(i,1:4) = [node_numbering, x, y, z];
%
%   nodecoor_cell is a 1-by-P cell array. nodecoor_cell{i} represents the 
%   node numbering and node coordinates in the i-th phase.
%       nodecoor_cell{i}(j,1:4) = [node_numbering, x, y, z];
%
%   ele_cell is a 1-by-P cell array. ele_cell{i} represent elements in the 
%   i-th phase. ele_cell{i} is an numeric array. Each row in ele_cell{i} is
%   the element numbering and node indices of an element. 
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    nodecoor_list = zeros( size(vert,1), 4 );
    nodecoor_list( :, 1 ) = 1: size(vert,1);
    nodecoor_list( :, 2:4 ) = vert;
    
    % extract the mesh of each phase from mesh
    [ nodecoor_cell, ele_cell ] = whole2phase( vert,ele,tnum );
    
end

function [ nodecoor_cell, ele_cell ] = whole2phase( vert, ele, tnum )
% whole2phase: extract the mesh of each phase from meshes
%
% input:
%   vert(k,1:3) = [x,y,z] coordinatex of k-th node 
%   ele(m,:) = [node_numbering_of_all_nodes] of m-th element
%   tnum(m,1) = n; means the m-th element is belong to phase n
%
% output:
%   nodecoor_cell{i} and ele_cell{i} represent phase i
%   nodecoor_cell{i}(j,1:4) = [node_numbering, x, y, z]
%   ele_cell{i}(j,:) = [element numbering, node_numbering_of_all_nodes]
%
    
    % phase label
    label_vec = unique( tnum );
    num_phase = length( label_vec );
    
    nodecoor_cell = cell( 1, num_phase );
    ele_cell = cell( 1, num_phase );

    for i = 1: num_phase
        %-----------------------------------------------------------------
        % get mesh element (node indices) of phase_i from ele
        tria_pI = ele( tnum==label_vec(i), : );
        
        % get unique node numbering of phase_i
        nodes_pI = unique( tria_pI(:) );
        
        %-----------------------------------------------------------------
        % extract node coordinates from vert according to node numbering
        % nodecoor_pI = [ node_numbering, x, y, z ]

        nodecoor_pI = zeros( length(nodes_pI), 4 );
        nodecoor_pI(:,1) = nodes_pI;       % node numbering
        
        for j = 1: size(nodecoor_pI,1)
              nodecoor_pI(j,2:4) = vert( nodecoor_pI(j,1), : );
        end

        %-----------------------------------------------------------------
        % get element numbering and node numbering of phase_i
        % ele_pI = [ element_numbering, node_numbering_in_a_element ]

        ele_pI = zeros( size(tria_pI,1), 1+size(ele,2) );
        ele_pI(:,1) = find( tnum==label_vec(i) );   % element numbering
        ele_pI(:,2:end) = tria_pI;                  % node numbering
        
        %-----------------------------------------------------------------
        nodecoor_cell{i} = nodecoor_pI;
        ele_cell{i} = ele_pI;
        
        %-----------------------------------------------------------------
    end
end










