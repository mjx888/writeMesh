function printInp3d( vert, ele, tnum, ele_type, precision, file_name, opt )
% printInp3d: write 3d finite element mesh (nodes and elements) to inp 
% file (Abaqus). Tested in software Abaqus 6.14. 
%
% Works for 3d mesh with single or multiple phases. Note that phase is also 
% known as part, domain, or physical surface. In finite element modeling of
% composite materials, each phase in the mesh represents a distinct 
% material component.
%
% Works for linear and quadratic element.
% Works for tetrahedral and hexahedral mesh.
%
% See the link below for usage examples of function printInp3d.
%   https://github.com/mjx888/writeMesh/blob/main/README.md
%
% Use functions: getNodeEle3d.m
%                getBCNode3d.m   getInterf3d.m
%
% usage:
%   printInp3d( vert, ele );
%   printInp3d( vert, ele, [], [], [], file_name );
%   printInp3d( vert, ele, tnum );
%   printInp3d( vert, ele, tnum, [], precision );
%   printInp3d( vert, ele, tnum, ele_type, precision );
%   printInp3d( vert, ele, tnum, ele_type, precision, file_name );
%   printInp3d( vert, ele, tnum, ele_type, precision, file_name, opt );
%
% input:
%   tnum, ele_type, precision, file_name, opt are optional.
%
%   vert: Mesh nodes. It’s a Nn-by-2 matrix, where 
%         Nn is the number of nodes in the mesh. Each row of vert 
%         contains the x, y coordinates for that mesh node.
%     
%   ele: Mesh elements. 
%        For linear tetrahedral elements, it s a Ne-by-4 matrix. 
%        For linear hexahedral elements, it s a Ne-by-8 matrix
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
%   opt - a structure array. It is the extra options for printInp3d.
%         It stores extra parameter settings for printInp3d.
%
%   opt.tf_printMaxMinNode - Boolean. Value: 0 or 1. Whether to print nodes
%                            at max & min location as node set.
%                            Default value: 1
%
%   opt.tf_printInterfNode - Boolean. Value: 0 or 1. Whether to print nodes
%                            at the interface as node set.
%                            Default value: 1
%
%   opt.tolerance - Numeric value. It's the tolerance for coordinates when
%                   searching boundary node set (at max & min location).
%                   Sub-routine will automatically find extrema for x y z 
%                   coordinates in the mesh. Those surface nodes with 
%                   coordinate satisfy |coordinate - extrema|< tolerance
%                   will be considered as boundary node set (at max & min 
%                   location). See writeMesh demo09 for usage example.
%                   Default value: 1E-10
%
%   opt.user_nodeSet - User-defined node set. A nested cell array. 
%                      Default value: {}
%             Exampe:  opt.user_nodeSet{1} = { 'name1', [2 5 8] };
%                      opt.user_nodeSet{2} = { 'name2', [36 23 56 80] };
%
%   opt.mode - Mode of printInp3d. This parameter is used to configure the 
%              arrangement of text within the inp file, or its format.
%              Value: 1, 2, or 3. Default value: 1
%     When =1, the declaration of assembly and instance will be neglected. 
%             This is the concise format.
%     When =2, assembly and instance will be declared explicitly (normal 
%             mode). The exported inp file would have a model with one 
%             part, which contains multiple sections. Each section 
%             corresponds to one phase in the mesh. 
%     When =3, assembly and instance will be declared explicitly (normal 
%             mode). The exported inp file would have a model with multiple
%             parts, where each part corresponds to one phase in the mesh.
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
        ele_wid = size( ele, 2 );
        
        switch ele_wid
            case 4
                ele_type = 'C3D4';
            case 10
                ele_type = 'C3D10';
            case 8
                ele_type = 'C3D8';
            case 20
                ele_type = 'C3D20';
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
    % fix node ordering for elements with negative volume
    % ---------------------------------------------------------------------
    % ele = fixOrdering3d( vert, ele );     % Under development
    
    % ---------------------------------------------------------------------
    % check mode & print inp
    % ---------------------------------------------------------------------
    
    if opt.mode == 1
        % Concise mode. Neglect the declaration of Assembly & Instance.
        % Create a model with one part, which contains multiple sections.
        printMode1( vert, ele, tnum, ele_type, precision, file_name, opt );
    
    elseif opt.mode == 2
        % Normal mode. Explicitly declare of assembly and instance.
        % Create a model with one part, which contains multiple sections.
        printMode2( vert, ele, tnum, ele_type, precision, file_name, opt );
    
    elseif opt.mode == 3
        % Normal mode. Explicitly declare of assembly and instance.
        % Create a model with multiple parts, where each part corresponds 
        % to one phase in the mesh.
        printMode3( vert, ele, tnum, ele_type, precision, file_name, opt );
        
    else
        error('Input "opt.mode" must be 1, 2, or 3.');
    end
    
    % ---------------------------------------------------------------------
end


function new_opt = setOption( opt )
% setOption: verify field names in opt and set values in new_opt according
% to opt

    % initialize new_opt with default field names & value 
    new_opt.tf_printMaxMinNode = true;
    new_opt.tf_printInterfNode = true;
    new_opt.tolerance = 1E-10;
    new_opt.user_nodeSet = {};
    new_opt.mode = 1;
    
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


function printMode1( vert, ele, tnum, ele_type, precision, file_name, opt )
% printMode1: Mode 1.
% Concise mode. Neglect the declaration of Assembly & Instance.
% Create a model with one part, which contains multiple sections.
    
	% format of inp file (if opt.mode=1)
	% ---------------------------------------------------------------------
	% Heading
	%
	% Node
    % Element
    % Section
    %
    % Node set
    %
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % prepare for writing file
    % ---------------------------------------------------------------------
    % Add node numbering and element numbering, and organize elements into 
    % cell array. eleC{i} represent elements in the i-th phase.
    
    [ nodecoor, nodecoorC, eleC ] = getNodeEle3d( vert, ele, tnum );
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
    % 3,4.69,23.82,0.50
    
    fprintf( fid, ...
        [ ...
        fmNodeNum, ',', fmNodeCo, ',', fmNodeCo, ',', fmNodeCo, '\n' ...
        ], ...
        nodecoor' ...
        );

    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    % Element
    
    for i = 1: num_phase
        % example:
        % *Element, type=CPS3, elset=Set-A

        fprintf( fid, [...
            '*Element, type=%s, elset=Set-%c'  '\n'...
            ], ele_type, num2char(i) );
        
        % example:
        % 3,173,400,475

        printEle( fid, eleC{i}, fmEleNum, fmNodeNum );

        fprintf( fid, '%s\n', '**' );
    end
    
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
        printNsMaxMin( fid, nodecoor, nodecoorC, ele, opt.tolerance );
    end
    
    % node set at the interface
    if opt.tf_printInterfNode
        printNsInterf( fid, nodecoorC );
    end
    
    % print user-defined node set
    if ~isempty( opt.user_nodeSet )
        printUserNs( fid, opt.user_nodeSet );
    end

    % ---------------------------------------------------------------------
    % u can use node set to define boundary condition
    % example:
    % ** Name: BC-FixY Type: Displacement/Rotation
    % *Boundary
    % Set-A-Ymin, 2, 2
    
    % ---------------------------------------------------------------------
    fclose(fid);
	
	disp('printInp3d Done! Check the inp file!');
    % ---------------------------------------------------------------------
end


function printMode2( vert, ele, tnum, ele_type, precision, file_name, opt )
% printMode2: Mode 2.
% Normal mode. Explicitly declare of assembly and instance.
% Create a model with one part, which contains multiple sections.

	% format of inp file (if opt.mode=2)
	% ---------------------------------------------------------------------
	% Heading
	%
    % Part
	%   Node
    %   Element
    %   Section
    %
    % Assembly
    %   Instance
    %
    %   Node set
    %
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % prepare for writing file
    % ---------------------------------------------------------------------
    % Add node numbering and element numbering, and organize elements into 
    % cell array. eleC{i} represent elements in the i-th phase.
    
    [ nodecoor, nodecoorC, eleC ] = getNodeEle3d( vert, ele, tnum );
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

    fprintf( fid, '%s\n', '**' );

	% ---------------------------------------------------------------------
    fprintf( fid, '%s\n', '*Part, name=Part-1' );
    
	% ---------------------------------------------------------------------
    % Node
    fprintf( fid, '*Node\n' );
    
    % print coordinates of nodes
    % example:
    % 3,4.69,23.82,0.50
    
    fprintf( fid, ...
        [ ...
        fmNodeNum, ',', fmNodeCo, ',', fmNodeCo, ',', fmNodeCo, '\n' ...
        ], ...
        nodecoor' ...
        );
    
    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    % Element
    
    for i = 1: num_phase
        % example:
        % *Element, type=CPS3, elset=Set-A

        fprintf( fid, [...
            '*Element, type=%s, elset=Set-%c'  '\n'...
            ], ele_type, num2char(i) );
        
        % example:
        % 3,173,400,475

        printEle( fid, eleC{i}, fmEleNum, fmNodeNum );

        fprintf( fid, '%s\n', '**' );
    end
    
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
    fprintf( fid, '%s\n', '*End Part' );
    fprintf( fid, '%s\n', '**' );
    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    fprintf( fid, '%s\n', '** ASSEMBLY' );
    fprintf( fid, '%s\n', '**' );
    fprintf( fid, '%s\n', '*Assembly, name=Assembly-1' );
    fprintf( fid, '%s\n', '**' );
    fprintf( fid, '%s\n', '*Instance, name=Instance-1, part=Part-1' );
    fprintf( fid, '%s\n', '*End Instance' );
    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    % print node set
    instanceName = 'Instance-1';

    % node set at max & min location
    if opt.tf_printMaxMinNode
        printNsMaxMin( fid, nodecoor, nodecoorC, ele, opt.tolerance, instanceName );
    end
    
    % node set at the interface
    if opt.tf_printInterfNode
        printNsInterf( fid, nodecoorC, instanceName );
    end
    
    % print user-defined node set
    if ~isempty( opt.user_nodeSet )
        printUserNs( fid, opt.user_nodeSet, instanceName );
    end

    % ---------------------------------------------------------------------
    fprintf( fid, '%s\n', '**' );
    fprintf( fid, '%s\n', '*End Assembly' );
    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    fclose(fid);
	
	disp('printInp3d Done! Check the inp file!');
    % ---------------------------------------------------------------------
end


function printMode3( vert, ele, tnum, ele_type, precision, file_name, opt )
% printMode3. Mode 3.
% Normal mode. Explicitly declare of assembly and instance.
% Create a model with multiple parts, where each part corresponds to one 
% phase in the mesh.

	% format of inp file (if opt.mode=3)
	% ---------------------------------------------------------------------
	% Heading
	%
    % Part-A
	%   Node
    %   Element
    %   Section
    %
    % Part-B
	%   Node
    %   Element
    %   Section
    %
    % Assembly
    %   Instance-A
    %   Instance-B
    %
    %   Node set
    %
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % prepare for writing file
    % ---------------------------------------------------------------------
    % Add node numbering and element numbering, and organize elements into 
    % cell array. eleC{i} represent elements in the i-th phase.
    
    [ ~, nodecoorC, eleC ] = getNodeEle3d( vert, ele, tnum );
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

    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    for i = 1: num_phase
        
	    % -----------------------------------------------------------------
        fprintf( fid, ['*Part, name=Part-%c' '\n'], num2char(i) );
        
	    % -----------------------------------------------------------------
        % Node
        fprintf( fid, '*Node\n' );
        
        % print coordinates of nodes
        % example:
        % 3,4.69,23.82,0.50
        
        fprintf( fid, ...
            [ ...
            fmNodeNum, ',', fmNodeCo, ',', fmNodeCo, ',', fmNodeCo, '\n' ...
            ], ...
            nodecoorC{i}' ...
            );
        
        fprintf( fid, '%s\n', '**' );
    
        % -----------------------------------------------------------------
        % Element
        
        % example:
        % *Element, type=CPS3, elset=Set-A

        fprintf( fid, [...
            '*Element, type=%s, elset=Set-%c'  '\n'...
            ], ele_type, num2char(i) );
        
        % example:
        % 3,173,400,475     % linear tria element

        printEle( fid, eleC{i}, fmEleNum, fmNodeNum );

        fprintf( fid, '%s\n', '**' );
        
        % -----------------------------------------------------------------
        % Section
    
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
        
        fprintf( fid, '%s\n', '**' );
    
	    % -----------------------------------------------------------------
        fprintf( fid, '%s\n', '*End Part' );
        fprintf( fid, '%s\n', '**' );
    end

    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    fprintf( fid, '%s\n', '** ASSEMBLY' );
    fprintf( fid, '%s\n', '**' );
    fprintf( fid, '%s\n', '*Assembly, name=Assembly-1' );
    fprintf( fid, '%s\n', '**' );

    for i = 1: num_phase
        fprintf( fid, ['*Instance, name=Instance-%c, part=Part-%c', '\n'], ...
                num2char(i), num2char(i) );

        fprintf( fid, '%s\n', '*End Instance' );
        fprintf( fid, '%s\n', '**' );
    end

    % ---------------------------------------------------------------------
    % print node set

    % node set at max & min location
    if opt.tf_printMaxMinNode
        printNsMaxMinXParts( fid, nodecoorC, ele, opt.tolerance );
    end
    
    % node set at the interface
    if opt.tf_printInterfNode
        printNsInterfXParts( fid, nodecoorC );
    end

    % ---------------------------------------------------------------------
    fprintf( fid, '%s\n', '**' );
    fprintf( fid, '%s\n', '*End Assembly' );
    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    fclose(fid);
	
	disp('printInp3d Done! Check the inp file!');
    % ---------------------------------------------------------------------
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
    % 3,173,400,475

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

function printNsMaxMin( fid, nodecoor, nodecoorC, ele, tolerance, instanceName )
% printNsMaxMin: print node set at max min location
% Used by printMode1, printMode2
%
% usage:
%   printNsMaxMin( fid, nodecoor, nodecoorC, ele, tolerance );               % for printMode1
%   printNsMaxMin( fid, nodecoor, nodecoorC, ele, tolerance, instanceName ); % for printMode2

    % ---------------------------------------------------------------------
    if nargin < 6
        instanceName = [];
    end
    
    num_phase = length(nodecoorC);
    
    % ---------------------------------------------------------------------
    % remove internal nodes in nodecoor, nodecoorC 
    
    % check
    if size(nodecoor, 1) > intmax('int32')
        error('Contact the author.')
    end
    
    surface_node = int32( getSurfaceNode3d(ele) );  % index
    
    % nodecoor
    pristine_node = int32( nodecoor(:,1) );  % index
    tf_vec = ismember( pristine_node, surface_node);
    nodecoor = nodecoor( tf_vec, : );

    % nodecoorC
    for i = 1: num_phase
        pristine_node = int32( nodecoorC{i}(:,1) );  % index
        tf_vec = ismember( pristine_node, surface_node);
        nodecoorC{i} = nodecoorC{i}( tf_vec, : );
    end
    
    % ---------------------------------------------------------------------
    % get node set at max min location

    % node set at xmin, xmax, ymin, ymax (globally)
    [ xmin_node, xmax_node, ...
      ymin_node, ymax_node, ...
      zmin_node, zmax_node ] = getBCNode3d( {nodecoor}, tolerance );
    
    xmin_node = xmin_node{1};
    xmax_node = xmax_node{1};
    ymin_node = ymin_node{1};
    ymax_node = ymax_node{1};
    zmin_node = zmin_node{1};
    zmax_node = zmax_node{1};

    % node set at xmin, xmax, ymin, ymax for each phase
    [ xmin_node_cell, xmax_node_cell, ...
      ymin_node_cell, ymax_node_cell, ...
      zmin_node_cell, zmax_node_cell ] = getBCNode3d( nodecoorC, tolerance );

    % ---------------------------------------------------------------------
    % node set at xmin, xmax, ymin, ymax, zmin, zmax (globally)
    % ---------------------------------------------------------------------
    % xmin
    if ~isempty( xmin_node )
        if isempty( instanceName )
            fprintf( fid, ['*Nset, nset=Set-Xmin' '\n'] );
        else
            fprintf( fid, ['*Nset, nset=Set-Xmin, instance=%s' '\n'], instanceName );
        end
        
	    printSet( fid, xmin_node );
    end

    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % xmax
    if ~isempty( xmax_node )
        if isempty( instanceName )
            fprintf( fid, ['*Nset, nset=Set-Xmax' '\n'] );
        else
            fprintf( fid, ['*Nset, nset=Set-Xmax, instance=%s' '\n'], instanceName );
        end

	    printSet( fid, xmax_node );
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % ymin
    if ~isempty( ymin_node )
        if isempty( instanceName )
            fprintf( fid, ['*Nset, nset=Set-Ymin' '\n'] );
        else
            fprintf( fid, ['*Nset, nset=Set-Ymin, instance=%s' '\n'], instanceName );
        end

	    printSet( fid, ymin_node );
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % ymax
    if ~isempty( ymax_node )
        if isempty( instanceName )
            fprintf( fid, ['*Nset, nset=Set-Ymax' '\n'] );
        else
            fprintf( fid, ['*Nset, nset=Set-Ymax, instance=%s' '\n'], instanceName );
        end
        
	    printSet( fid, ymax_node );
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % zmin
    if ~isempty( zmin_node )
        if isempty( instanceName )
            fprintf( fid, ['*Nset, nset=Set-Zmin' '\n'] );
        else
            fprintf( fid, ['*Nset, nset=Set-Zmin, instance=%s' '\n'], instanceName );
        end

	    printSet( fid, zmin_node );
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % zmax
    if ~isempty( zmax_node )
        if isempty( instanceName )
            fprintf( fid, ['*Nset, nset=Set-Zmax' '\n'] );
        else
            fprintf( fid, ['*Nset, nset=Set-Zmax, instance=%s' '\n'], instanceName );
        end
        
	    printSet( fid, zmax_node );
    end
    
    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    % node set at xmin, xmax, ymin, ymax, zmin, zmax for each phase
    % ---------------------------------------------------------------------
    if num_phase == 1
        return
    end
    % ---------------------------------------------------------------------
    % xmin
    for i = 1: num_phase
	    if ~isempty( xmin_node_cell{i} )
            if isempty( instanceName )
                fprintf( fid, ['*Nset, nset=Set-Xmin-%c' '\n'], num2char(i) );
            else
                fprintf( fid, ['*Nset, nset=Set-Xmin-%c, instance=%s' '\n'], ...
                        num2char(i), instanceName );
            end
    
		    printSet( fid, xmin_node_cell{i} );
	    end
    end

    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % xmax
    for i = 1: num_phase
	    if ~isempty( xmax_node_cell{i} )
            if isempty( instanceName )
                fprintf( fid, ['*Nset, nset=Set-Xmax-%c' '\n'], num2char(i) );
            else
                fprintf( fid, ['*Nset, nset=Set-Xmax-%c, instance=%s' '\n'], ...
                        num2char(i), instanceName );
            end
    
		    printSet( fid, xmax_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % ymin
    for i = 1: num_phase
	    if ~isempty( ymin_node_cell{i} )
            if isempty( instanceName )
                fprintf( fid, ['*Nset, nset=Set-Ymin-%c' '\n'], num2char(i) );
            else
                fprintf( fid, ['*Nset, nset=Set-Ymin-%c, instance=%s' '\n'], ...
                        num2char(i), instanceName );
            end
    
		    printSet( fid, ymin_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % ymax
    for i = 1: num_phase
	    if ~isempty( ymax_node_cell{i} )
            if isempty( instanceName )
                fprintf( fid, ['*Nset, nset=Set-Ymax-%c' '\n'], num2char(i) );
            else
                fprintf( fid, ['*Nset, nset=Set-Ymax-%c, instance=%s' '\n'], ...
                        num2char(i), instanceName );
            end
    
		    printSet( fid, ymax_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % zmin
    for i = 1: num_phase
	    if ~isempty( zmin_node_cell{i} )
            if isempty( instanceName )
                fprintf( fid, ['*Nset, nset=Set-Zmin-%c' '\n'], num2char(i) );
            else
                fprintf( fid, ['*Nset, nset=Set-Zmin-%c, instance=%s' '\n'], ...
                        num2char(i), instanceName );
            end
    
		    printSet( fid, zmin_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % zmax
    for i = 1: num_phase
	    if ~isempty( zmax_node_cell{i} )
            if isempty( instanceName )
                fprintf( fid, ['*Nset, nset=Set-Zmax-%c' '\n'], num2char(i) );
            else
                fprintf( fid, ['*Nset, nset=Set-Zmax-%c, instance=%s' '\n'], ...
                        num2char(i), instanceName );
            end
    
		    printSet( fid, zmax_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
end


function printNsInterf( fid, nodecoorC, instanceName )
% printNsInterf: print node set at the interface
% Used by printMode1, printMode2
%
% usage:
%   printNsInterf( fid, nodecoorC );                % for printMode1
%   printNsInterf( fid, nodecoorC, instanceName );  % for printMode2

    % ---------------------------------------------------------------------
    if nargin < 3
        instanceName = [];
    end

    % ---------------------------------------------------------------------
    num_phase = length(nodecoorC);
    interfnode_cell = getInterf3d( nodecoorC );
    % interfnode_cell{i,j} are nodes at interface i,j

    for i = 1: num_phase-1
	    for j = i+1: num_phase
		    if ~isempty( interfnode_cell{i,j} )
                if isempty( instanceName )
                    fprintf( fid, ['*Nset, nset=Set-Interf-%c%c' '\n'], ...
                            num2char(i), num2char(j) );
                else
                    fprintf( fid, ['*Nset, nset=Set-Interf-%c%c, instance=%s' '\n'], ...
                            num2char(i), num2char(j), instanceName );
                end
    
			    printSet( fid, interfnode_cell{i,j} );
		    end
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
end

function printUserNs( fid, nodeSet, instanceName )
% printUserNs: print user-defined node set
%
% usage:
%   printUserNs( fid, nodeSet )                 % for printMode1
%   printUserNs( fid, nodeSet, instanceName );  % for printMode2

    % ---------------------------------------------------------------------
    if nargin < 3
        instanceName = [];
    end
    
    % ---------------------------------------------------------------------
    num_set = length(nodeSet);
    
    for i = 1: num_set
	    if ~isempty( nodeSet{i}{1} ) && ~isempty( nodeSet{i}{2} )
            
            if isempty( instanceName )
                fprintf( fid, ['*Nset, nseprintNsMaxMint=Set-%s' '\n'], nodeSet{i}{1} );
            else
                fprintf( fid, ['*Nset, nset=Set-%s, instance=%s' '\n'], ...
                        nodeSet{i}{1}, instanceName );
            end
    
		    printSet( fid, nodeSet{i}{2} );

            fprintf( fid, '%s\n', '**' );
	    end
    end
    % ---------------------------------------------------------------------
end

function printNsMaxMinXParts( fid, nodecoorC, ele, tolerance )
% printNsMaxMinXParts: print node set at max min location for each part
% Used by function printMode3
% 

    % ---------------------------------------------------------------------
    % remove internal nodes in nodecoorC 
    
    num_phase = length(nodecoorC);
    surface_node = int32( getSurfaceNode3d(ele) );  % index
    
    for i = 1: num_phase
        % check
        if size(nodecoorC{i}, 1) > intmax('int32')
            error('Contact the author.')
        end

        pristine_node = int32( nodecoorC{i}(:,1) );  % index
        tf_vec = ismember( pristine_node, surface_node);
        nodecoorC{i} = nodecoorC{i}( tf_vec, : );
    end

    % ---------------------------------------------------------------------
    % get node set

    % node set at xmin, xmax, ymin, ymax for each phase
    [ xmin_node_cell, xmax_node_cell, ...
      ymin_node_cell, ymax_node_cell, ...
      zmin_node_cell, zmax_node_cell ] = getBCNode3d( nodecoorC, tolerance );

    % ---------------------------------------------------------------------
    % node set at xmin, xmax, ymin, ymax, zmin, zmax for each phase
    % ---------------------------------------------------------------------
%     if num_phase == 1
%         return
%     end
    
    % ---------------------------------------------------------------------
    % xmin
    for i = 1: num_phase
	    if ~isempty( xmin_node_cell{i} )
            fprintf( fid, ...
                    ['*Nset, nset=Set-Xmin-%c, instance=Instance-%c' '\n'], ...
                    num2char(i), num2char(i) );
            
		    printSet( fid, xmin_node_cell{i} );
	    end
    end

    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % xmax
    for i = 1: num_phase
	    if ~isempty( xmax_node_cell{i} )
            fprintf( fid, ...
                    ['*Nset, nset=Set-Xmax-%c, instance=Instance-%c' '\n'], ...
                    num2char(i), num2char(i) );

		    printSet( fid, xmax_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );

    % ---------------------------------------------------------------------
    % ymin
    for i = 1: num_phase
	    if ~isempty( ymin_node_cell{i} )
            fprintf( fid, ...
                    ['*Nset, nset=Set-Ymin-%c, instance=Instance-%c' '\n'], ...
                    num2char(i), num2char(i) );
    
		    printSet( fid, ymin_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % ymax
    for i = 1: num_phase
	    if ~isempty( ymax_node_cell{i} )
            fprintf( fid, ...
                    ['*Nset, nset=Set-Ymax-%c, instance=Instance-%c' '\n'], ...
                    num2char(i), num2char(i) );
    
		    printSet( fid, ymax_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    
    % ---------------------------------------------------------------------
    % zmin
    for i = 1: num_phase
	    if ~isempty( zmin_node_cell{i} )
            fprintf( fid, ...
                    ['*Nset, nset=Set-Zmin-%c, instance=Instance-%c' '\n'], ...
                    num2char(i), num2char(i) );
    
		    printSet( fid, zmin_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
    % zmax
    for i = 1: num_phase
	    if ~isempty( zmax_node_cell{i} )
            fprintf( fid, ...
                    ['*Nset, nset=Set-Zmax-%c, instance=Instance-%c' '\n'], ...
                    num2char(i), num2char(i) );
    
		    printSet( fid, zmax_node_cell{i} );
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
end


function printNsInterfXParts( fid, nodecoorC )
% printNsInterfXParts: print node set at the interface for each part
% Used by function printMode3
%

    % ---------------------------------------------------------------------
    num_phase = length(nodecoorC);
    interfnode_cell = getInterf3d( nodecoorC );
    % interfnode_cell{i,j} are nodes in part i at interface i,j
    % interfnode_cell{j,i} are nodes in part j at interface i,j
    
    for i = 1: num_phase-1
	    for j = i+1: num_phase
		    if ~isempty( interfnode_cell{i,j} )
                % nodes in part i at interface i,j
                fprintf( fid, ...
                    ['*Nset, nset=Set-Interf-%c%c-in-%c, instance=Instance-%c' '\n'], ...
                        num2char(i), num2char(j), num2char(i), num2char(i) );
                
			    printSet( fid, interfnode_cell{i,j} );
                
                % nodes in part j at interface i,j
                fprintf( fid, ...
                    ['*Nset, nset=Set-Interf-%c%c-in-%c, instance=Instance-%c' '\n'], ...
                        num2char(i), num2char(j), num2char(j), num2char(j) );
                
			    printSet( fid, interfnode_cell{j,i} );
		    end
	    end
    end
    
    fprintf( fid, '%s\n', '**' );
    % ---------------------------------------------------------------------
end













