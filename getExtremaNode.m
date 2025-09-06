function node_cell = getExtremaNode( nodecoor_cell, xyzchar, fHandle, tolerance )
% getExtremaNode: get nodes with extrema coordinates, such as x min, x max
% xyzchar - 'x', 'y', or 'z'
%
% example:
%     xyzchar = 'z';
%     tolerance = 0.1;
%     f = @min;
%     node_cell = getExtrema( nodecoor_cell, xyzchar, f, tolerance );
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % --------------------------------------------------------------------
    switch xyzchar
        case 'x'
            idx = 2;
        case 'y'
            idx = 3;
        case 'z'
            idx = 4;
        otherwise
            error('not correct input - xyzchar')
    end

    % --------------------------------------------------------------------
    % find global extrema
    num_phase = length( nodecoor_cell );
    extrema_vec = zeros( 1, num_phase );
    for i = 1: num_phase
        extrema_vec(i) = fHandle( nodecoor_cell{i}(:,idx) );
    end
    extrema = fHandle( extrema_vec );
    
    % --------------------------------------------------------------------
    % find node index for extrema
    node_cell = cell( 1, num_phase );
    
    for i = 1: num_phase
        nc = nodecoor_cell{i};
        tf_vec = ( abs( nc(:,idx)-extrema ) < tolerance );
        node_cell{i} = nc( tf_vec, 1 );  % record node index
    end
    % --------------------------------------------------------------------
end






















