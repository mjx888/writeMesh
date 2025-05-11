function printInp2d( vert, ele, tnum, ele_type, precision, file_name, opt )
% printInp2d: write 2d finite element mesh (nodes and elements) to inp 
%           file (Abaqus). Tested in software Abaqus. 
%           The exported inp file will have a model with one part, which 
%           contains multiple sections. Each section corresponds to one 
%           material phase in the mesh.
%           
%           Works for linear and quadratic element.
%           Works for triangular and quadrilateral element.
%           Function printInp2d will automatically export node set.
%           opt.user_nodeSet is used to defined customized node set.
%
%           Use functions: getNodeEle.m  fixOrdering.m  
%                          getBCNode.m   getInterf.m
%
% usage:
%   printInp2d( vert, ele );
%   printInp2d( vert, ele, [], [], [], file_name );
%   printInp2d( vert, ele, tnum );
%   printInp2d( vert, ele, tnum, [], precision );
%   printInp2d( vert, ele, tnum, ele_type, precision );
%   printInp2d( vert, ele, tnum, ele_type, precision, file_name );
%   printInp2d( vert, ele, tnum, ele_type, precision, file_name, opt );
%
% input:
%   tnum, ele_type, precision, file_name, opt are optional.
%
%   vert: Mesh nodes. It’s a Nn-by-2 matrix, where 
%         Nn is the number of nodes in the mesh. Each row of vert 
%         contains the x, y coordinates for that mesh node.
%     
%   ele: Mesh elements. 
%        For linear triangular elements, it s a Ne-by-3 matrix. 
%        For linear quadrilateral elements, it s a Ne-by-4 matrix

%        Ne is the number of elements in the mesh. Each row in ele 
%        contains the indices of the nodes for that mesh element.
%   
%   tnum: Label of phase, which corresponds to physical surface tag in Gmsh. 
%         tnum is a Ne-by-1 array, where Ne is the number of elements.
%         tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%         When omitted, assign one phase.
%     
%   ele_type: element type in Abaqus.
%             When omitted, ele_typ will be determined based on the size of 
%             variable ele.
%
%   precision: number of digits to the right of the decimal point 
%              when writing node coordinates.
%              When omitted, precision=8;
%
%   file_name: file name of inp file, such as 'aaa.inp', 'D:\aaa.inp'.
%              When omitted, file_name='test.inp';
%
%   opt - a structure array. It is the extra options for printInp2d.
%         It stores extra parameter settings for printInp2d.
%
%   opt.tf_printMaxMinNode - Boolean. Value: 0 or 1. Whether to print nodes
%                            at max & min location as node set.
%                            Default value: 1
%
%   opt.tf_printInterfNode - Boolean. Value: 0 or 1. Whether to print nodes
%                            at the interface as node set.
%                            Default value: 1
%
%   opt.user_nodeSet - User-defined node set. A nested cell array. 
%                      Default value: {}
%             Exampe:  opt.user_nodeSet{1} = { 'name1', [2 5 8] };
%                      opt.user_nodeSet{2} = { 'name2', [36 23 56 80] };
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
%                  https://github.com/mjx888/writeMesh
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
    % Node set
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
        precision = [];
    end
    
    if nargin < 6
        file_name = [];
    end

    if nargin < 7
        opt = [];
    end

    % ---------------------------------------------------------------------
    % verify field names and set values for opt
    opt = setOption( opt );

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
    
    if isempty(precision)
        precision = 8;
    end

    if isempty(file_name)
        % write to current folder
        file_name = 'test.inp';
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
    
    [ nodecoor, nodecoorC, eleC ] = getNodeEle( vert, ele, tnum );
    
    num_phase = length( eleC );
    
    % ---------------------------------------------------------------------
    % format of number

    % format_node_coor
    % '%.(precision)f'
    fmNodeCo = [ '%.', num2str( precision ), 'f' ];    
    
    fmNodeNum = '%d';     % format_node_num
    fmEleNum = '%d';      % format_ele_num
    
    % ---------------------------------------------------------------------
    % start writing to file
    % ---------------------------------------------------------------------
	fid=fopen(file_name,'wW');
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
    
    for i = 1: num_phase
        % example:
        % *Element, type=CPS3, elset=Set-A
        fprintf( fid, [...
            '*Element, type=%s, elset=Set-%c'  '\n'...
            ], ele_type, num2char(i) );
        
        % example:
        % 3,173,400,475     % linear tria element
        % 87,428,584,561,866,867,868    % quadratic tria element

        printEle( fid, eleC{i}, fmEleNum, fmNodeNum );
    end
    
    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    % Section

    for i = 1: num_phase
        % example:
        % ** Section: Section-A
        % *Solid Section, elset=Set-A, material=Material-A
        % ,

        fprintf( fid, [...
            '** Section: Section-%c'            '\n'...
            '*Solid Section, elset=Set-%c, material=Material-%c'  '\n'...
            ','                                 '\n'...
            ], ...
            num2char(i), ...
            num2char(i), num2char(i) );
    end
    
    fprintf( fid, '%s\n', '**' );
    
    % ---------------------------------------------------------------------
    % print node set

    % node set at max & min location
    if opt.tf_printMaxMinNode
        printNodeSetMaxMin( fid, nodecoor, nodecoorC );
    end
    
    % node set at the interface
    if opt.tf_printInterfNode
        printNodeSetInterface( fid, nodecoorC );
    end
    
    % print user-defined node set
    if ~isempty( opt.user_nodeSet )
        printUserNodeSet( fid, opt.user_nodeSet );
    end

    % ---------------------------------------------------------------------
    % u can use node set to define boundary condition
    % example:
    % ** Name: BC-FixY Type: Displacement/Rotation
    % *Boundary
    % Set-A-Ymin, 2, 2

    % ---------------------------------------------------------------------
    fclose(fid);
	
	disp('printInp2d Done! Check the inp file!');
    % ---------------------------------------------------------------------
end


function new_opt = setOption( opt )
% setOption: verify field names in opt and set values in new_opt according
% to opt

    % initialize new_opt with default field names & value 
    new_opt.tf_printMaxMinNode = true;
    new_opt.tf_printInterfNode = true;
    new_opt.user_nodeSet = {};
    
    if isempty(opt)
        return
    end

    if ~isstruct(opt)
        error("opt is not a structure array. Not valid input.")
    end

    % get the field names of opt
    nameC = fieldnames(opt);

    % verify field names in opt and set values in new_opt
    % compare the field name of opt with new_opt using for loop
    % if a field name of opt exist in new_opt, assign the that field value 
    % in opt to new_opt
    % if a field name of opt not exist in new_opt, show error

    for i = 1: length(nameC)
        if isfield( new_opt, nameC{i} )
            value = getfield( opt, nameC{i} );
            new_opt = setfield( new_opt, nameC{i}, value );
        else
            error("Field name %s in opt is not correct.", nameC{i});
        end
    end

end


function phase_char = num2char( k )
% num2char: convert number 1 2 3 to character A B C

    phase_char = char( k-1+65 );     % 'ABCD...'
end

function printEle( fid, ele, format_ele_num, format_node_num )
% printEle: print elements
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

function printSet( fid, nodeSet )
% printSet: print node set

    for i=1:length(nodeSet)
        if mod( i, 16 ) == 0 || i == length(nodeSet)
            fprintf( fid, '%d\n', nodeSet(i) );
        else
            fprintf( fid, '%d, ', nodeSet(i) );
        end
    end
end

function printNodeSetMaxMin( fid, nodecoor, nodecoorC )
% printNodeSetMaxMin: print node set at max min location

    % ---------------------------------------------------------------------
    % get node set

    num_phase = length(nodecoorC);

    [ xmin_node_cell, xmax_node_cell, ...
      ymin_node_cell, ymax_node_cell ] = getBCNode( nodecoorC );
    
    [ xmin_node, xmax_node, ...
      ymin_node, ymax_node ] = getBCNode( {nodecoor} );
    xmin_node = xmin_node{1};
    xmax_node = xmax_node{1};
    ymin_node = ymin_node{1};
    ymax_node = ymax_node{1};
    
    % ---------------------------------------------------------------------
    % node set at xmin, xmax, ymin, ymax (globally)
    % ---------------------------------------------------------------------
    % xmin
    if ~isempty( xmin_node )
	    fprintf( fid, [...
		    '*Nset, nset=Set-Xmin'   '\n'...
		    ] );

	    printSet( fid, xmin_node );
    end

    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % xmax
    if ~isempty( xmax_node )
	    fprintf( fid, [...
		    '*Nset, nset=Set-Xmax'   '\n'...
		    ] );

	    printSet( fid, xmax_node );
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % ymin
    if ~isempty( ymin_node )
	    fprintf( fid, [...
		    '*Nset, nset=Set-Ymin'   '\n'...
		    ] );

	    printSet( fid, ymin_node );
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % ymax
    if ~isempty( ymax_node )
	    fprintf( fid, [...
		    '*Nset, nset=Set-Ymax'   '\n'...
		    ] );
        
	    printSet( fid, ymax_node );
    end
    
    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    % node set at xmin, xmax, ymin, ymax for each phase
    % ---------------------------------------------------------------------
    % xmin
    for i = 1: num_phase
	    if ~isempty( xmin_node_cell{i} )
		    fprintf( fid, [...
			    '*Nset, nset=Set-Xmin-%c'   '\n'...
			    ], num2char(i) );
    
		    printSet( fid, xmin_node_cell{i} );
	    end
    end

    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % xmax
    for i = 1: num_phase
	    if ~isempty( xmax_node_cell{i} )
		    fprintf( fid, [...
			    '*Nset, nset=Set-Xmax-%c'   '\n'...
			    ], num2char(i) );
    
		    printSet( fid, xmax_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % ymin
    for i = 1: num_phase
	    if ~isempty( ymin_node_cell{i} )
		    fprintf( fid, [...
			    '*Nset, nset=Set-Ymin-%c'   '\n'...
			    ], num2char(i) );
    
		    printSet( fid, ymin_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % ymax
    for i = 1: num_phase
	    if ~isempty( ymax_node_cell{i} )
		    fprintf( fid, [...
			    '*Nset, nset=Set-Ymax-%c'   '\n'...
			    ], num2char(i) );
    
		    printSet( fid, ymax_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
end


function printNodeSetInterface( fid, nodecoorC )
% printNodeSetInterface: print node set at the interface

    num_phase = length(nodecoorC);
    interfnode_cell = getInterf( nodecoorC );

    for i = 1: num_phase-1
	    for j = i+1: num_phase
		    if ~isempty( interfnode_cell{i,j} )
			    % node at interface i,j
			    fprintf( fid, [...
				    '*Nset, nset=Set-Interf-%c%c' '\n'...
				    ], num2char(i), num2char(j) );
    
			    printSet( fid, interfnode_cell{i,j} );
		    end
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
end

function printUserNodeSet( fid, nodeSet )
% printUserNodeSet: print user-defined node set

    num_set = length(nodeSet);
    
    for i = 1: num_set
	    if ~isempty( nodeSet{i}{1} ) && ~isempty( nodeSet{i}{2} )
		    fprintf( fid, [...
			    '*Nset, nset=Set-%s'   '\n'...
			    ], nodeSet{i}{1} );
    
		    printSet( fid, nodeSet{i}{2} );

            fprintf( fid, '%s\n', '**' );
	    end
    end
    
end




