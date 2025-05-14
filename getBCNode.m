function [ xmin_node_cell, xmax_node_cell, ...
           ymin_node_cell, ymax_node_cell ] = getBCNode( nodecoor_cell )
% getBCNode: get node set related to boundary condition (BC)
% e.g., find xmin first, then get node number, whose x coordinate is equal 
%       to xmin
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    num_phase = length( nodecoor_cell );
    
    % ---------------------------------------------------------------------
    % xmin_node_cell
    xmin_vector = zeros( 1, num_phase );  % xmin_vector stores xmin for each phase
    for i = 1: num_phase
        xmin_vector(i) = min( nodecoor_cell{i}(:,2) );
    end
    xmin = min( xmin_vector );      % global xmin
    idx = find( xmin_vector == xmin );  % means xmin happened in which phase
    
    xmin_node_cell = cell( 1, num_phase );
    for i = idx
        xmin_node_cell{i} = getXYNode( nodecoor_cell{i}, xmin, 'x' );
    end

    % ---------------------------------------------------------------------
    % xmax_node_cell
    xmax_vector = zeros( 1, num_phase );
    for i = 1: num_phase
        xmax_vector(i) = max( nodecoor_cell{i}(:,2) );
    end
    xmax = max( xmax_vector );
    idx = find( xmax_vector == xmax );
    
    xmax_node_cell = cell( 1, num_phase );
    for i = idx
        xmax_node_cell{i} = getXYNode( nodecoor_cell{i}, xmax, 'x' );
    end
    
    % ---------------------------------------------------------------------
    % ymin_node_cell
    ymin_vector = zeros( 1, num_phase );
    for i = 1: num_phase
        ymin_vector(i) = min( nodecoor_cell{i}(:,3) );
    end
    ymin = min( ymin_vector );
    idx = find( ymin_vector == ymin );
    
    ymin_node_cell = cell( 1, num_phase );
    for i = idx
        ymin_node_cell{i} = getXYNode( nodecoor_cell{i}, ymin, 'y' );
    end
    
    % ---------------------------------------------------------------------
    % ymax_node_cell
    ymax_vector = zeros( 1, num_phase );
    for i = 1: num_phase
        ymax_vector(i) = max( nodecoor_cell{i}(:,3) );
    end
    ymax = max( ymax_vector );
    idx = find( ymax_vector == ymax );
    
    ymax_node_cell = cell( 1, num_phase );
    for i = idx
        ymax_node_cell{i} = getXYNode( nodecoor_cell{i}, ymax, 'y' );
    end
    % ---------------------------------------------------------------------
end


function xy_node_phase = getXYNode( nodecoor_phase, value, xychar )
% find node in nodecoor_phase
% return node number, whose x or y coordinate is equal to value
% xychar - 'x' or 'y'

    switch xychar
        case 'x'
            idx = 2;
        case 'y'
            idx = 3;
        otherwise
            error('not correct input')
    end
    
    thresh = 1E-10;  % threshold
    xy_node_phase = zeros( length(nodecoor_phase), 1 );
    
    for i = 1: length(nodecoor_phase)
        if abs( nodecoor_phase(i,idx)-value ) < thresh
            xy_node_phase(i) = nodecoor_phase( i, 1 );
        end
    end
    xy_node_phase( xy_node_phase==0 ) = [];
    
end
