function plotMeshes( vert, ele, tnum )
% plotMeshes: plot triangular finite element mesh
% Also works for quadratic or quadrilateral elements
%
% Nodes must be counter-clockwise ordering in an linear element.
%
% usage:
%   plotMeshes( vert, ele, tnum );  % multiple phases
%   plotMeshes( vert, ele );        % one phase
%
% input:
%   verrt - Node data. N-by-2 array.
%       vert(i,1:2) = [x_coordinate, y_coordinate] of the i-th node
%
%   ele - Node numbering for each element. 
%       For example, if linear triangle element, ele is M-by-3 array.
%       ele(j,1:3) = [node_numbering_of_3_nodes] of the j-th element
%
%   tnum - Label of material phase. P-by-1 array.
%       tnum(j,1) = k; means the j-th element is belong to the k-th phase
%
%
% Im2mesh is copyright (C) 2019-2025 by Jiexian Ma and is distributed under
% the terms of the GNU General Public License (version 3).
% 
% Project website: https://github.com/mjx888/im2mesh
%

    %--------------------------------------------------------------------
    % check the number of inputs
    if nargin == 2
        tnum = ones(size(ele,1),1);
    elseif nargin == 3
        % normal case
    else
        error("check the number of inputs");
    end
    
    %--------------------------------------------------------------------
    % check element type
    ele_wid = size(ele,2);

    if ele_wid == 3         % linear triangle
        range_vec = 1:3;
    elseif ele_wid == 6     % quadratic triangle
        range_vec = [1 4 2 5 3 6];
    elseif ele_wid == 4     % linear quadrilateral
        range_vec = 1:4;
    elseif ele_wid == 8     % quadratic quadrilateral
        range_vec = [1 5 2 6 3 7 4 8];
    else
        error("ele - wrong size")
    end

    %--------------------------------------------------------------------
    % plot mesh
    figure;
    hold on; 
    axis image off;

    tvalue = unique( tnum );
    num_phase = length( tvalue );
    
    % setup color
    if num_phase == 1
        col = 0.98;
    elseif num_phase > 1
        col = 0.3: 0.68/(num_phase-1): 0.98;
    else
        error("num_phase < 1")
    end
    
    % use function patch to plot
    for i = 1: num_phase
        current_phase = tvalue(i);
        patch( ...
            'faces',ele( tnum==current_phase, range_vec ), ...
            'vertices',vert, ...
            'facecolor',[ col(i), col(i), col(i) ], ...
            'edgecolor',[.1,.1,.1] ...
            );
    end
    hold off
    %--------------------------------------------------------------------
end

