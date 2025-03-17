function [ nodecoor_list, nodecoor_cell, ele_cell ] = getNodeEle( vert, tria, tnum, ~ )
% getNodeEle: get node coordinares and elements from mesh
%
% Detail: Add node numbering and element numbering, and organize elements 
% into cell array. eleC{i} represent elements in the i-th phase.
%
% Works for linear and quadratic element.
% Works for triangular and quadrilateral element.
%
% usage:
%   [ nodecoor_list, nodecoor_cell, ele_cell ] = getNodeEle( vert, tria, tnum );
%
% input:
%   verrt - Node data. N-by-2 array.
%       vert(i,1:2) = [x_coordinate, y_coordinate] of the i-th node
%
%   tria - Node numbering for each triangle. 
%       M-by-3 array (linear element) or M-by-6 array (2nd order element)
%       tria(j,1:3) = [node_numbering_of_3_nodes] of the j-th element
%       tria(j,1:6) = [node_numbering_of_6_nodes] of the j-th element
%
%   tnum - Label of material phase. P-by-1 array.
%       tnum(j,1) = k; means the j-th element is belong to the k-th phase
%
% output:
%   nodecoor_list is node coordinates, N-by-3 array.
%       nodecoor_list(i,1:3) = [node_numbering, x_coordinate, y_coordinate];
%
%   nodecoor_cell is a 1-by-P cell array. nodecoor_cell{i} represents the 
%   node numbering and node coordinates in the i-th phase.
%       nodecoor_cell{i}(j,1:3) = [node_numbering, x_coordinate, y_coordinate];
%
%   ele_cell is a 1-by-P cell array. ele_cell{i} represent elements in the 
%   i-th phase.
%       % When the elements are linear, ele_cell{i} has 4 rows.
%       ele_cell{i}(j,1:4) = [ element_numbering, node_numbering_of_node_1, ...
%                                               node_numbering_of_node_2, ...
%                                               node_numbering_of_node_3 
%                             ];
%
%       % When the elements are second order, ele_cell{i} has 7 rows.
%       ele_cell{i}(j,1:7) = [ element_numbering, node_numbering_of_all nodes ]
%
%
% Copyright (C) 2019-2025 by Jiexian Ma
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % check the number of inputs
    if nargin == 3
        % normal case
    elseif nargin == 4
        warning("the 4th parameter - ele_order is no longer support.")
    else
        error('check the number of inputs');
    end
    
    nodecoor_list = zeros( size(vert,1), 3 );
    nodecoor_list( :, 1 ) = 1: size(vert,1);
    nodecoor_list( :, 2:3 ) = vert;
    
    % extract the mesh of each phase from mesh
    [ nodecoor_cell, ele_cell ] = whole2phase( vert,tria,tnum );
    
end

function [ nodecoor_cell, ele_cell ] = whole2phase( vert, tria, tnum )
% whole2phase: extract the mesh of each phase from meshes
%
% input:
%   vert(k,1:2) = [x_coordinate, y_coordinate] of k-th node 
%   tria(m,:) = [node_numbering_of_all_nodes] of m-th element (3 or 6 nodes)
%   tnum(m,1) = n; means the m-th element is belong to phase n
%
% output:
%   nodecoor_cell{i} and ele_cell{i} represent phase i
%   nodecoor_cell{i}(j,1:3) = [node_numbering, x_coordinate, y_coordinate]
%   ele_cell{i}(j,1:4/7) = [element numbering, node_numbering_of_all_nodes]
%
    
    % phase label
    label_vec = unique( tnum );
    num_phase = length( label_vec );
    
    nodecoor_cell = cell( 1, num_phase );
    ele_cell = cell( 1, num_phase );

    for i = 1: num_phase
        %-----------------------------------------------------------------
        % get triangular mesh (node numbering of 3 nodes) of phase_i from tria
        tria_pI = tria( tnum==label_vec(i), : );
        
        % get unique node numbering of phase_i
        nodes_pI = unique( tria_pI(:) );
        
        %-----------------------------------------------------------------
        % extract node coordinates from vert according to node numbering
        % nodecoor_pI = [ node_numbering, x, y ]

        nodecoor_pI = zeros( length(nodes_pI), 3 );
        nodecoor_pI(:,1) = nodes_pI;       % node numbering
        
        for j = 1: size(nodecoor_pI,1)
              nodecoor_pI(j,2:3) = vert( nodecoor_pI(j,1), : );
        end

        %-----------------------------------------------------------------
        % get element numbering and node numbering of phase_i
        % ele_pI = [ element_numbering, node_numbering_in_a_element ]

        ele_pI = zeros( size(tria_pI,1), 1+size(tria,2) );
        ele_pI(:,1) = find( tnum==label_vec(i) );   % element numbering
        ele_pI(:,2:end) = tria_pI;                  % node numbering
        
        %-----------------------------------------------------------------
        nodecoor_cell{i} = nodecoor_pI;
        ele_cell{i} = ele_pI;

    end
end










