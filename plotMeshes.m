function plotMeshes( vert, ele, tnum, color_code )
% plotMeshes: plot triangular or quadrilateral finite element mesh
% Works for linear and quadratic elements
%
% Nodes must be counter-clockwise ordering in an linear element.
%
% usage:
%   plotMeshes( vert, ele );        % one phase
%   plotMeshes( vert, ele, tnum );  % multiple phases
%
%   plotMeshes( vert, ele, [], color_code );    % one phase
%   plotMeshes( vert, ele, tnum, color_code );  % multiple phases
%   plotMeshes( vert, ele, tnum, 2 );
%
% input:
%   vert - Node data. N-by-2 array.
%       vert(i,1:2) = [x_coordinate, y_coordinate] of the i-th node
%
%   ele - Node numbering for each element. 
%       For example, if linear triangle element, ele is M-by-3 array.
%       ele(j,1:3) = [node_numbering_of_3_nodes] of the j-th element
%
%   tnum - Label of material phase. P-by-1 array.
%       tnum(j,1) = k; means the j-th element is belong to the k-th phase
%
%   color_code - Color code for selecting colormap.
%                Interger. Value: 0-10. Default value: 0.
%
% by Jiexian Ma, mjx0799@gmail.com
% 
% Project website: https://github.com/mjx888/im2mesh
%

    %--------------------------------------------------------------------
    % check the number of inputs
    if nargin == 2
        tnum = ones(size(ele,1),1);
        color_code = 0;

    elseif nargin == 3
        if isempty(tnum)
            tnum = ones(size(ele,1),1);
        end
        color_code = 0;

    elseif nargin == 4
       if isempty(tnum)
            tnum = ones(size(ele,1),1);
       end

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
    tvalue = unique( tnum );
    num_phase = length( tvalue );
    
    %--------------------------------------------------------------------
    % setup color
    % Create variable 'colors' - num_phase-by-3 array.
    % Each row in 'colors' is one rgb color.

    switch color_code
        case 0
            % grayscale
            if num_phase == 1
                col = 0.98;
            elseif num_phase > 1
                col = 0.3: 0.68/(num_phase-1): 0.98;
                col = col(:);
            end
            colors = [col, col, col];

        case 1
            colors = lines( num_phase );
        case 2
            colors = parula( num_phase );
        case 3
            colors = turbo( num_phase );
        case 4
            colors = jet( num_phase );
        case 5
            colors = hot( num_phase );
        case 6
            colors = cool( num_phase );
        case 7
            colors = summer( num_phase );
        case 8
            colors = winter( num_phase );
        case 9
            colors = bone( num_phase );
        case 10
            colors = pink( num_phase );
        otherwise
            error('Input argument color_code is out of range.')
    end

    %--------------------------------------------------------------------
    % plot mesh
    figure;
    hold on; 
    axis image off;
    
    % use function patch to plot
    for i = 1: num_phase
        current_phase = tvalue(i);
        patch( ...
            'faces', ele( tnum==current_phase, range_vec ), ...
            'vertices', vert, ...
            'facecolor', colors(i,:), ...
            'edgecolor', [.1,.1,.1] ...
            );
    end
    hold off
    
    %--------------------------------------------------------------------
end

