function ele = fixOrdering( vert, ele )
% fixOrdering: fix node ordering in each elements of 2D finite 
% element mesh. Node ordering in a element should be counterclockwise.
%
% Detail: Caculate area sign of an element. If negative, flip ordering of
%         nodes in a element. Counterclockwise ordering produce positive 
%         area. Clockwise ordering produce negative area.
%
% usage:
%   ele = fixOrdering( vert, ele );
% 
% input:
%   vert: Mesh nodes. It’s a Nn-by-2 matrix, where 
%         Nn is the number of nodes in the mesh. Each row of vert 
%         contains the x, y coordinates for that mesh node.
%     
%   ele: Mesh elements. For linear triangular elements, 
%         it s a Ne-by-3 matrix, where Ne is the number of elements in 
%         the mesh. Each row in ele contains the indices of the nodes 
%         for that mesh element.
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% 
% Project website: https://github.com/mjx888/im2mesh
%
    
    vert = vert(:,1:2);
    ele_width = size(ele,2);

    switch ele_width
        case 3
            area = triarea( vert, ele );
            ele(area<0.,:) = ele( area<0., [1,3,2] );
        case 6
            area = triarea( vert, ele );
            ele(area<0.,:) = ele( area<0., [1,3,2,6,5,4] );
        case 4
            area = quadarea( vert, ele );
            ele(area<0.,:) = ele( area<0., [1,4,3,2] );
        case 8
            area = quadarea( vert, ele );
            ele(area<0.,:) = ele( area<0., [1,4,3,2,8,7,6,5] );
        otherwise
            error("fixOrdering only work for triangle and quadrangle.");
    end

    numNegative = sum( area<0. );
    if numNegative > 0
        disp([' ', num2str(numNegative), ' negative area elements are fixed.'])
    end
end

function [area] = triarea(vert,tria)
% triarea: calculates the signed area of each triangular element

    ev12 = vert(tria(:,2),:)-vert(tria(:,1),:) ;
    ev13 = vert(tria(:,3),:)-vert(tria(:,1),:) ;

    area = ev12(:,1).*ev13(:,2) ...
         - ev12(:,2).*ev13(:,1) ;
    area = 0.5 * area;
end

function area = quadarea(vert, ele)
% quadarea: calculates the signed area of each quadrilateral element,

    % Number of elements
    Ne = size(ele, 1);
    
    area = zeros(Ne, 1);

    for i = 1:Ne
        % Extract the node indices for element i
        nodes = ele(i, :);
        
        % Extract the coordinates of the 4 vertices
        x = vert(nodes, 1);
        y = vert(nodes, 2);
        
        % Shoelace formula components
        sum1 = x(1)*y(2) + x(2)*y(3) + x(3)*y(4) + x(4)*y(1);
        sum2 = y(1)*x(2) + y(2)*x(3) + y(3)*x(4) + y(4)*x(1);

        % Signed area: no absolute value
        area(i) = 0.5 * (sum1 - sum2);
    end
end

