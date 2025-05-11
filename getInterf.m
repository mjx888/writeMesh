function interfnode_cell = getInterf( nodecoor_cell )
% get interface node between phase nodecoor_cell{i} and nodecoor_cell{j}
%
% interfnode_cell - n*n cell array. If nodecoor_cell{i} and
%                   nodecoor_cell{j} don't contact,  interfnode_cell{i,j}
%                   would be []. interfnode_cell{i,j} stores interfacial
%                   node number of nodecoor_cell{i} .
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% 
% Project website: https://github.com/mjx888/im2mesh
%

    len = length( nodecoor_cell );
    interfnode_cell = cell( len, len );
    
    for i = 1: len-1
        for j = i+1: len
            nodecoor_A = nodecoor_cell{i};
            nodecoor_B = nodecoor_cell{j};
            [interfnode_A, interfnode_B] = findInterfNode(nodecoor_A, nodecoor_B);
            
            if isempty(interfnode_A) || isempty(interfnode_B)
                continue
            else
                interfnode_cell{i,j} = interfnode_A;
                interfnode_cell{j,i} = interfnode_B;
            end
        end
    end

end

function [interfnode_A, interfnode_B] = findInterfNode(nodecoor_A, nodecoor_B)
% find the interface between phase A and phase B, return node number in 
% corresponding phase
    
    mark_A = false( length(nodecoor_A), 1 );
    mark_B = false( length(nodecoor_B), 1 );
    
    for i = 1: size( nodecoor_A( :, 2:3 ), 1 )
        idx = find( nodecoor_A( i, 2 ) == nodecoor_B( :, 2 ) & ...
                    nodecoor_A( i, 3 ) == nodecoor_B( :, 3 ) );
        len = length(idx);

        switch len
            case 0
                continue
            case 1
                mark_A( i ) = 1;
                mark_B( idx ) =1;
            otherwise
                disp('multiple contact points')
        end
    end

    interfnode_A = nodecoor_A( mark_A, 1 ); % find 1
    interfnode_B = nodecoor_B( mark_B, 1 );

end











