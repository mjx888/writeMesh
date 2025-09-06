function printBdf3d( vert, ele, tnum, ele_type, precision, file_name )
% printBdf3d: write 3d finite element mesh (nodes and elements) to bdf 
% file (Nastran bulk data, compatible with COMSOL). 
%           
% Works for linear tetrahedral and linear hexahedral element.
% Not work for quadratic element.
%
% Use functions: getNodeEle3d.m
%
% See the link below for usage examples.
%   https://github.com/mjx888/writeMesh/blob/main/README.md
%
% usage:
%   printBdf3d( vert, ele );
%   printBdf3d( vert, ele, [], [], [], file_name );
%   printBdf3d( vert, ele, tnum );
%   printBdf3d( vert, ele, tnum, [], precision );
%   printBdf3d( vert, ele, tnum, [], precision, file_name );
%
% input:
%   ele_type, precision, file_name are optional.
%
%   vert: Mesh nodes. It’s a Nn-by-3 matrix, where 
%         Nn is the number of nodes in the mesh. Each row of vert 
%         contains the x, y, z coordinates for that mesh node.
%     
%   ele: Mesh elements. 
%        For linear tetrahedral elements, it s a Ne-by-4 matrix. 
%        For linear hexahedral elements, it s a Ne-by-8 matrix
%         
%        Ne is the number of elements in the mesh. Each row in ele 
%        contains the indices of the nodes for that mesh element.
%   
%   tnum: Label of phase, which corresponds to physical surface tag in Gmsh. 
%         tnum is a Ne-by-1 array, where Ne is the number of elements.
%         tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%         When omitted, assign one phase.
%     
%   ele_type: Please set ele_type as an empty array.
%             ele_type is reserved for future development.
%
%   precision: number of digits to the right of the decimal point 
%              when writing node coordinates.
%              When omitted, precision=8;
%
%   file_name: file name of bdf file, such as 'aaa.bdf', 'D:\aaa.bdf'.
%              When omitted, file_name='test.bdf';
%
%
% This is sub-project of Im2mesh package. If you use this function, please
% cite as follows: 
%  Ma, J., & Li, Y. (2025). Im2mesh: A MATLAB/Octave package for generating
%  finite element mesh based on 2D multi-phase image (2.1.5). Zenodo. 
%  https://doi.org/10.5281/zenodo.14847059
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
% 
% Project website: https://github.com/mjx888/im2mesh
%                  https://github.com/mjx888/writeMesh
%

    % format of bdf file
    % ---------------------------------------------------------------------
    % BEGIN BULK
    % GRID*,1,,0.50000000,0.50000000,*
    % *,0.00000000
    % GRID*,2,,0.50000000,3.50000000,*
    % *,0.00000000
    % CTETRA*,1,1,628,208,*
    % *,59,246
    % CTETRA*,2,1,665,671,*
    % *,285,178
    % ENDDATA
    
    % ---------------------------------------------------------------------
    % check inputs
    % ---------------------------------------------------------------------
    % Check the number of inputs. If missing, set as empty. 
    if nargin < 2
        error("Not enough input arguments.");
    end
    
    if nargin < 3
        tnum = [];
    end
    
    if nargin < 4
        ele_type = [];
    end
    
    if nargin < 5
        precision = [];
    end
    
    if nargin < 6
        file_name = [];
    end

    % ---------------------------------------------------------------------
    % check input size
    if size(vert,2) >= 4
        warning("Only 1st to 3rd column in vert will be considered.");
        vert = vert( :, 1:3 );
    end
    
    if ~isempty(tnum) && size(tnum,1) ~= size(ele,1)
        error("The 3rd input argument tnum has wrong size.");
    end
    % ---------------------------------------------------------------------
    % If input is empty, assign defaualt value to input
    if isempty(tnum)
        tnum = ones( size(ele,1), 1 );
    end
    
    if isempty(ele_type)
        % do nothing
    end
    
    if isempty(precision)
        precision = 8;
    end

    if isempty(file_name)
        % write to current folder
        file_name = 'test.bdf';
    end
    
    % ---------------------------------------------------------------------
    % check input type

    % Validate that 'precision' is a positive integer
    a = precision;
    if ~isnumeric(a) || ~isscalar(a) || a <= 0 || mod(a, 1) ~= 0
        error('Input "precision" must be a positive integer.');
    end

    % Validate that 'file_name' is a string
    b = file_name;
    if ~(ischar(b) || isstring(b))
        error('Input "file_name" must be a string.');
    end

    % ---------------------------------------------------------------------
    % fix node ordering for elements with negative area
    % ele = fixOrdering3d( vert, ele );     % Under development
    
    % ---------------------------------------------------------------------
    % prepare for writing file
    % ---------------------------------------------------------------------
    % Add node numbering and element numbering, and organize elements into 
    % cell array. eleC{i} represent elements in the i-th phase.

    [nodecoor, ~, eleC] = getNodeEle3d( vert, ele, tnum );

    % ---------------------------------------------------------------------
    numNode = size( nodecoor, 1 );

    % ---------------------------------------------------------------------
    % format of number
    
    % field width of node numbering
    width_node_num = 1 + floor( log10( numNode ) );         % 18964 -> 5
    if width_node_num > 16
        error('more than 16 digits')
    end
    
    % num_digits_of_int_part
    numDigitsIntPart = 1 + floor( log10( max(nodecoor(end,2:4)) ) );
                                                             % 182.9 -> 3
    if (numDigitsIntPart + precision + 1) > 16
        error('more than 16 digits')
    end
    
    % format_node_coor
    % '%.(precision)f'
    fmNodeCo = [ '%.', num2str( precision ), 'f' ];     

    % ---------------------------------------------------------------------
    % start writing to file
    % ---------------------------------------------------------------------
	fid=fopen( file_name, 'wW' );
    % ---------------------------------------------------------------------
    fprintf( fid, 'BEGIN BULK\n');

    fprintf( fid, [...
        '$ BDF file generated by Im2mesh package'              '\n'...
        ]);
    
    % print node
    % GRID*,3,,0.5000,2.5000,*
    % *,3.5000

    fprintf( fid, ...
            [ 'GRID*,%d,,', fmNodeCo, ',', fmNodeCo, ',*\n', ...
             '*,', fmNodeCo, '\n'], ...
            nodecoor' ...
            );

    % ---------------------------------------------------------------------
    % renumber global element numbering in eleC{i}(:,1)
    num_phase = length( eleC );
    count = 0;
    for i = 1: num_phase
        eleC{i}(:,1) = (1:size(eleC{i},1))' + count;
        count = count + size(eleC{i},1);
    end

    % ---------------------------------------------------------------------
    % print element
    ele_wid = size( eleC{1}, 2 ) -1;
    
    if ele_wid == 4
        % linear tetrahedral element
        % CTETRA*,5,1,40,46,*
        % *,47,73

        for i = 1: num_phase
            fprintf( fid, ...
                ['CTETRA*,%d,%d,%d,%d', ',*\n', ...
                 '*,', '%d,%d\n'], ...
                [ eleC{i}(:,1), i * ones(size(eleC{i},1),1), eleC{i}(:,2:5) ]' ...
                );
        end
        
    elseif ele_wid == 8
        % linear hexahedral element
        % CHEXA*,5,1,40,46,*
        % *,47,41,10,16,*
        % *,17,11
        
        for i = 1: num_phase
            fprintf( fid, ...
                ['CHEXA*,%d,%d,%d,%d', ',*\n', ...
                 '*,', '%d,%d,%d,%d', ',*\n', ...
                 '*,', '%d,%d\n'], ...
                [ eleC{i}(:,1), i * ones(size(eleC{i},1),1), eleC{i}(:,2:9) ]' ...
                );
        end
    else
        error('Function printBdf3d do not support quadratic elements.');
    end
    
    % ---------------------------------------------------------------------
	fprintf( fid, 'ENDDATA' );
    
    % ---------------------------------------------------------------------
    fclose(fid);
	
	disp('printBdf3d Done! Check the bdf file!');
    % ---------------------------------------------------------------------
end
