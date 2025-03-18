function printInp2d( vert, ele, tnum, ele_type, precision_nodecoor, path_file_name )
% printInp2d: write 2d finite element mesh (nodes and elements) to inp 
%           file (Abaqus). Test in software Abaqus. 
%           The exported inp file will have a model with one part, which 
%           contains multiple sections. Each section corresponds to one 
%           material phase in the mesh.
%           Use functions: getNodeEle.m  fixOrdering.m
%
%           Works for linear and quadratic element.
%           Works for triangular and quadrilateral element.
%
% usage:
%   printInp2d( vert, ele );
%   printInp2d( vert, ele, [], [], [], path_file_name );
%   printInp2d( vert, ele, tnum );
%   printInp2d( vert, ele, tnum, [], precision_nodecoor );
%   printInp2d( vert, ele, tnum, ele_type, precision_nodecoor );
%   printInp2d( vert, ele, tnum, ele_type, precision_nodecoor, path_file_name );
%
% input:
%   ele_type, precision_nodecoor, path_file_name are optional.
%
%   vert: Mesh nodes. It’s a Nn-by-2 matrix, where 
%         Nn is the number of nodes in the mesh. Each row of vert 
%         contains the x, y coordinates for that mesh node.
%     
%   ele: Mesh elements. 
%        For linear triangular elements, it s a Ne-by-3 matrix. 
%        For linear quadrilateral elements, it s a Ne-by-4 matrix
%         
%        Ne is the number of elements in the mesh. Each row in ele 
%        contains the indices of the nodes for that mesh element.
%   
%   tnum: Label of phase, which corresponds to physical surface tag in Gmsh. 
%         tnum is a Ne-by-1 array, where Ne is the number of elements.
%         tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%         When omitted, assign one phase.
%     
%   ele_type: element type in Abaqus
%             When omitted, ele_typ will be determined based on the size of 
%             variable ele.
%
%   precision_nodecoor: number of digits to the right of the decimal point 
%                       when writing node coordinates.
%                       When omitted, precision_nodecoor=8;
%
%   path_file_name: file name of inp file, such as 'aaa.inp', 'D:\aaa.inp'.
%                   When omitted, path_file_name='test.inp';
%
%
% This is sub-project of Im2mesh package. If you use this function, please
% cite as follows: 
%  Ma, J., & Li, Y. (2025). Im2mesh: A MATLAB/Octave package for generating
%  finite element mesh based on 2D multi-phase image (2.1.5). Zenodo. 
%  https://doi.org/10.5281/zenodo.14847059
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% 
% Project website: https://github.com/mjx888/im2mesh
%


	% format of inp file
	% ---------------------------------------------------------------------
	% Heading
	%
	% Node
    %
    % Element
    %
    % Section
    %
    % ---------------------------------------------------------------------

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
        precision_nodecoor = [];
    end
    
    if nargin < 6
        path_file_name = [];
    end

    % ---------------------------------------------------------------------
    % check input size
    if size(vert,2) >= 3
        warning("Z coordnates of mesh nodes will be ignored.");
        vert = vert( :, 1:2 );
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
        ele_wid = size( ele, 2 );
        
        switch ele_wid
            case 3
                ele_type = 'CPS3';
            case 6
                ele_type = 'CPS6';
            case 4
                ele_type = 'CPS4';
            case 8
                ele_type = 'CPS8';
            otherwise
                error("Bad input in mesh elements.")
        end
    end
    
    if isempty(precision_nodecoor)
        precision_nodecoor = 8;
    end

    if isempty(path_file_name)
        % write to current folder
        path_file_name = 'test.inp';
    end

    % ---------------------------------------------------------------------
    % check input type

    % Validate that 'precision_nodecoor' is a positive integer
    a = precision_nodecoor;
    if ~isnumeric(a) || ~isscalar(a) || a <= 0 || mod(a, 1) ~= 0
        error('Input "precision_nodecoor" must be a positive integer.');
    end

    % Validate that 'path_file_name' is a string
    b = path_file_name;
    if ~(ischar(b) || isstring(b))
        error('Input "path_file_name" must be a string.');
    end

    % Validate that 'ele_type' is a string
    c = ele_type;
    if ~(ischar(c) || isstring(c))
        error('Input "ele_type" must be a string.');
    end

    % ---------------------------------------------------------------------
    % fix node ordering for elements with negative area
    ele = fixOrdering( vert, ele );
    
    % ---------------------------------------------------------------------
    % prepare for writing file
    % ---------------------------------------------------------------------
    % Add node numbering and element numbering, and organize elements into 
    % cell array. eleC{i} represent elements in the i-th phase.
    
    [ nodecoor, ~, eleC ] = getNodeEle( vert, ele, tnum );
    
    % ---------------------------------------------------------------------
    % convert number 1 2 3 to character A B C
    num_sect = length( eleC );
    sect_ascii = 65: ( 65 + num_sect - 1);
    % section character
    sect_char = char( sect_ascii );     % 'ABCD...'
    
    % ---------------------------------------------------------------------
    % format of number

    % format_node_coor
    % '%.(precision)f'
    fmNodeCo = [ '%.', num2str( precision_nodecoor ), 'f' ];    
    
    fmNodeNum = '%d';     % format_node_num
    fmEleNum = '%d';      % format_ele_num
    
    % ---------------------------------------------------------------------
    % start writing to file
    % ---------------------------------------------------------------------
	fid=fopen(path_file_name,'wW');
    % ---------------------------------------------------------------------
	% Heading
    fprintf( fid, [...
        '*Heading'                                              '\n'...
        '*Preprint, echo=NO, model=NO, history=NO, contact=NO'  '\n'...
        '** INP file generated by Im2mesh package'              '\n'...
        '**'                                                    '\n'...
        ] ...
        );
    
	% ---------------------------------------------------------------------
    % Node
    fprintf( fid, '*Node\n' );
    
    % print coordinates of nodes
    % example:
    % 3,4.69000000,23.82000000
    %
    % '%d,%.4f,%.4f,%.4f\n'
    
    fprintf( fid, ...
        [ ...
        fmNodeNum, ',', fmNodeCo, ',', fmNodeCo, '\n' ...
        ], ...
        nodecoor' ...
        );
    
    % ---------------------------------------------------------------------
    % Element
    
    for i = 1: num_sect
        % example:
        % *Element, type=CPS3, elset=Set-A
        fprintf( fid, [...
            '*Element, type=%s, elset=Set-%c'  '\n'...
            ], ele_type, sect_char(i) );
        
        % example:
        % 3,173,400,475     % linear tria element
        % 87,428,584,561,866,867,868    % quadratic tria element

        printEle( fid, eleC{i}, fmEleNum, fmNodeNum );
    end
    
    % ---------------------------------------------------------------------
    % Section

    for i = 1: num_sect
        % example:
        % ** Section: Section-A
        % *Solid Section, elset=Set-A, material=Material-A
        % ,

        fprintf( fid, [...
            '** Section: Section-%c'            '\n'...
            '*Solid Section, elset=Set-%c, material=Material-%c'  '\n'...
            ','                                 '\n'...
            ], ...
            sect_char(i), ...
            sect_char(i), sect_char(i) );
    end
    
	fprintf( fid, '**' );
    
    % ---------------------------------------------------------------------
    fclose(fid);
	
	disp('printInp2d Done! Check the inp file!');
    % ---------------------------------------------------------------------
end

function printEle( fid, ele, format_ele_num, format_node_num )
% work for linear element and quadratic element

    num_column = size( ele, 2 );
    
    % fprintf( fid, '%d,%d,%d,%d,%d\n', ele' );
    % example:
    % 3,173,400,475                 % linear tria element
    % 87,428,584,561,866,867,868    % quadratic tria element

    fprintf( fid, ...
        [   format_ele_num, ',', ...
            repmat([format_node_num, ','], [1,num_column-2]), ...
            format_node_num, '\n' ...
        ], ...
        ele' ...
        );
end

% % old version
% function printEle( fid, ele, format_ele_num, format_node_num )
% % work for linear element and quadratic element
% 
%     num_column = size( ele, 2 );
% 
%     switch num_column
%         case 4
%             % linear element
%             % example:
%             % 3,173,400,475
%             fprintf( fid, ...
%                 [   format_ele_num, ',', ...
%                     repmat([format_node_num, ','], [1,2]), ...
%                     format_node_num, '\n' ...
%                 ], ...
%                 ele' ...
%                 );
%             
%         case 7
%             % quadratic element
%             % example:
%             % 87,428,584,561,866,867,868
%             fprintf( fid, ...
%                 [   format_ele_num, ',', ...
%                     repmat([format_node_num, ','], [1,5]), ...
%                     format_node_num, '\n' ...
%                 ], ...
%                 ele' ...
%                 );
%             
%         otherwise
%             disp('unidentified data')
%     end
% 
% end
