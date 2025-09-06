function node_set = displayBCNode3d( vert, ele, tnum, tolerance, markerSize, plane )
% displayBCNode3d: display boundary node set (at max & min location)
% See writeMesh demo09 for examples.
%
% usage:
%   displayBCNode3d( vert, ele, [], tolerance, markerSize, plane );
%   displayBCNode3d( vert, ele, tnum, tolerance, markerSize, plane );
%
% input:
%   tolerance - tolerance for coordinates when searching boundary node set 
%               (at max & min location). Sub-routine will automatically 
%               find extrema for x y z coordinates in the mesh. Those 
%               surface nodes with coordinate satisfy 
%               |coordinate - extrema|< tolerance will be considered as 
%               boundary node set (at max & min location). 
%
%   markerSize - marker size when plotting boundary node set
%
%   plane - value: 'all', 'x', 'y', or 'z'
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % ---------------------------------------------------------------------
    % If tnum is empty, assign defaualt value to tnum
    if isempty(tnum)
        tnum = ones( size(ele,1), 1 );
    end

    [ nodecoor, nodecoorC, ~ ] = getNodeEle3d( vert, ele, tnum );
    num_phase = length(nodecoorC);
    
    % ---------------------------------------------------------------------
    % remove internal nodes in nodecoorC 
    
    % check
    if size(nodecoor, 1) > intmax('int32')
        error('Contact the author.')
    end
    
    surface_node = int32( getSurfaceNode3d(ele) );  % index
    
    for i = 1: num_phase
        pristine_node = int32( nodecoorC{i}(:,1) );  % index
        tf_vec = ismember( pristine_node, surface_node);
        nodecoorC{i} = nodecoorC{i}( tf_vec, : );
    end

    %--------------------------------------------------------------------
    % find boundary node set

    [ xmin_node_cell, xmax_node_cell, ...
      ymin_node_cell, ymax_node_cell, ...
      zmin_node_cell, zmax_node_cell ] = getBCNode3d( nodecoorC, tolerance );
    
    switch plane
        case 'all'
            node_set = [
                        vertcat( xmin_node_cell{:} );
                        vertcat( xmax_node_cell{:} );
                        vertcat( ymin_node_cell{:} );
                        vertcat( ymax_node_cell{:} );
                        vertcat( zmin_node_cell{:} );
                        vertcat( zmax_node_cell{:} );
                        ];
        case 'x'
           node_set = [
                        vertcat( xmin_node_cell{:} );
                        vertcat( xmax_node_cell{:} );
                        ];
        case 'y'
           node_set = [
                        vertcat( ymin_node_cell{:} );
                        vertcat( ymax_node_cell{:} );
                        ];
        case 'z'
           node_set = [
                        vertcat( zmin_node_cell{:} );
                        vertcat( zmax_node_cell{:} );
                        ];
        otherwise
            error('not correct input - xyzchar')
    end
    
    %--------------------------------------------------------------------
    tvalue = unique( tnum );
    num_phase = length( tvalue );

    % set color
    if num_phase == 1
        col = 0.98;
    elseif num_phase > 1
        col = 0.3: 0.68/(num_phase-1): 0.98;
        col = col(:);
    end
    colors = [col, col, col];
    
    %--------------------------------------------------------------------
    % plot mesh
    figure;
    hold on;
    axis image off;
    set(gcf,'GraphicsSmoothing','off');
    
    % use function patch to plot faces and edges
    for i = 1: num_phase
        current_phase = tvalue(i);
        ele_temp = ele( tnum==current_phase, : );
        f_temp = getBoundaryFaces( ele_temp );

        patch( ...
            'faces', f_temp, ...
            'vertices', vert, ...
            'facecolor', colors(i,:), ...
            'edgecolor', [.1,.1,.1], ...
            'linewidth', 0.5, ...
            'facealpha', 1, ...
            'edgealpha', 0.5 ...
            );
    end
    
    % show node set
    for i = 1: length( node_set )
        idx = node_set(i);
        x = nodecoor( nodecoor(:,1)==idx, 2 );
        y = nodecoor( nodecoor(:,1)==idx, 3 );
        z = nodecoor( nodecoor(:,1)==idx, 4 );
        plot3( x,y,z, 'r.', 'MarkerSize', markerSize );
    end

    view([45 30])
    hold off
    %--------------------------------------------------------------------
end

function F = getBoundaryFaces(ele)
% getBoundaryFaces: find external facets of a mesh
%
% Inputs
%     ele - Ne×4 tetra connectivity
% Outputs
%     F - Nf×3 boundary triangle node indices 

    % ---------------------------------------------------------------------
    ele_wid = size( ele, 2 );
    
    if ele_wid == 10
        ele = ele( :, 1:4 );
    end
    if ele_wid == 20
        ele = ele( :, 1:8 );
    end
    
    % ---------------------------------------------------------------------
    ele_wid = size( ele, 2 );

    if ele_wid == 4
        % Four faces per tet
        faces  = [ ele(:,[2 3 4]);   
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

    % ---------------------------------------------------------------------
end


