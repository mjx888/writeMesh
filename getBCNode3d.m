function [ xmin_node_cell, xmax_node_cell, ...
           ymin_node_cell, ymax_node_cell, ...
           zmin_node_cell, zmax_node_cell ] = getBCNode3d( nodecoor_cell, tolerance )
% getBCNode3d: get node set for boundary condition (BC)
% 
% input:
%   tolerance - tolerance for coordinates when searching boundary node set 
%               (at max & min location). Sub-routine will automatically 
%               find extrema for x y z coordinates in the mesh. Those 
%               surface nodes with coordinate satisfy 
%               |coordinate - extrema|< tolerance will be considered as 
%               boundary node set (at max & min location). 
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % ---------------------------------------------------------------------
    % xmin_node_cell
    xyzchar = 'x';
    f = @min;
    xmin_node_cell = getExtremaNode( nodecoor_cell, xyzchar, f, tolerance );
    
    % ---------------------------------------------------------------------
    % xmax_node_cell
    xyzchar = 'x';
    f = @max;
    xmax_node_cell = getExtremaNode( nodecoor_cell, xyzchar, f, tolerance );
    
    % ---------------------------------------------------------------------
    % ymin_node_cell
    xyzchar = 'y';
    f = @min;
    ymin_node_cell = getExtremaNode( nodecoor_cell, xyzchar, f, tolerance );
    
    % ---------------------------------------------------------------------
    % ymax_node_cell
    xyzchar = 'y';
    f = @max;
    ymax_node_cell = getExtremaNode( nodecoor_cell, xyzchar, f, tolerance );
    
    % ---------------------------------------------------------------------
    % zmin_node_cell
    xyzchar = 'z';
    f = @min;
    zmin_node_cell = getExtremaNode( nodecoor_cell, xyzchar, f, tolerance );
    
    % ---------------------------------------------------------------------
    % zmax_node_cell
    xyzchar = 'z';
    f = @max;
    zmax_node_cell = getExtremaNode( nodecoor_cell, xyzchar, f, tolerance );
    
    % ---------------------------------------------------------------------
end









