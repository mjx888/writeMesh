function printMsh( vert, ele, tnum, conn, precision_nodecoor, path_file_name )
% printMsh: write 2d finite element mesh (nodes and elements) to msh file.
%           msh is Gmsh mesh file format. MSH file format version: 4.1
%           Test in software Gmsh 4.13.1
%           Use functions: tria2Surface.m  tria2BoundEdge.m  fixOrdering.m
%           
%           printMsh only works for 2d trangles & linear element.
%
%
% usage:
%   printMsh( vert, ele );
%   printMsh( vert, ele, [], [], [], path_file_name );
%   printMsh( vert, ele, tnum );
%   printMsh( vert, ele, [], [], precision_nodecoor );
%   printMsh( vert, ele, tnum, [], precision_nodecoor );
%   printMsh( vert, ele, tnum, [], precision_nodecoor );
%   printMsh( vert, ele, tnum, [], precision_nodecoor, path_file_name );
%
% input:
%   tnum, conn, precision_nodecoor, path_file_name are optional.
%
%   vert: Mesh nodes. It’s a Nn-by-2 matrix, where 
%         Nn is the number of nodes in the mesh. Each row of vert 
%         contains the x, y coordinates for that mesh node.
%     
%   ele: Mesh elements. For linear triangular elements, 
%         it s a Ne-by-3 matrix, where Ne is the number of elements in 
%         the mesh. Each row in ele contains the indices of the nodes 
%         for that mesh element.
%     
%   tnum: Label of phase, which corresponds to physical surface tag in Gmsh. 
%         tnum is a Ne-by-1 array, where Ne is the number of elements.
%         tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%         When omitted, assign one phase.
%     
%   conn: C-by-2 array of constraining edges. Each row defines an edge
%         If you don't know what is conn, just set conn as an empty array, 
%         which will not affect the generated msh file. When conn is an 
%         empty array, function printMsh will do computation to obatin the 
%         missing conn parameter.
%
%   precision_nodecoor: number of digits to the right of the decimal point 
%                       when writing node coordinates.
%                       When omitted, precision_nodecoor = 8;
%
%   path_file_name: file name of msh file, such as 'aaa.msh', 'D:\aaa.msh'.
%                   When omitted, path_file_name = 'test.msh';
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

    % format of msh file
    % ---------------------------------------------------------------------
    % $MeshFormat
    % 4.1 0 8
    % $EndMeshFormat
    % $Entities
    % 
    % $EndEntities
    % $Nodes
    % 
    % $EndNodes
    % $Elements
    % 
    % $EndElements
    
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
        conn = [];
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
    
    if size(ele,2) == 6
        warning("Mesh elements will be processed as linear trangles.");
        ele = ele( :, 1:3 );
    end
        
    if size(ele,2) ~= 3 && size(ele,2) ~= 6
        error("printMsh not work for quadrilateral elements.");
    end
    
    if ~isempty(tnum) && size(tnum,1) ~= size(ele,1)
        error("The 3rd input argument tnum has wrong size.");
    end
    
    % ---------------------------------------------------------------------
    % If input is empty, assign defaualt value to input
    if isempty(tnum)
        tnum = ones( size(ele,1), 1 );
    end
    
    if isempty(conn)
        conn = tria2BoundEdge( ele, tnum );
    end

    if isempty(precision_nodecoor)
        precision_nodecoor = 8;
    end

    if isempty(path_file_name)
        % write to current folder
        path_file_name = 'test.msh';
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

    % ---------------------------------------------------------------------
    % fix node ordering for elements with negative area
    ele = fixOrdering( vert, ele );
    
    % ---------------------------------------------------------------------
    % prepare for writing file
    % ---------------------------------------------------------------------
    % format of number
    % '%.(precision)f'
    
    % format_node_coor
    fmtNodeCo = [ '%.', num2str( precision_nodecoor ), 'f' ];
    
    fmtNodeNum = '%d';  % format_node_num
    fmtEleNum = '%d';   % format_ele_num

    % ---------------------------------------------------------------------
    % start writing to file
    % ---------------------------------------------------------------------
	fid=fopen(path_file_name,'wW');
    % ---------------------------------------------------------------------
	% Step1. Heading
    % ---------------------------------------------------------------------
    fprintf( fid, [...
        '$Comments'                                 '\n'...
        'MSH file generated by Im2mesh package'     '\n'...
        '$EndComments'                              '\n'...
        ] ...
        );

    % $MeshFormat
    % 4.1 0 8
    fprintf( fid, [...
        '$MeshFormat'       '\n'...
        '4.1 0 8'           '\n'...
        '$EndMeshFormat'    '\n'...
        ] ...
        );

    % ---------------------------------------------------------------------
    % Step2. Define Entities (point, curve, surface, physical surface)
    % ---------------------------------------------------------------------
    % get point, curve, surface for define physical surface
    % Here the points is for boundary representation. Not mesh node.
    
    % conn is the constraining edges
    % numPoints is the number of points on constraining edges
    pntIndx = unique(conn);
    pntIndx = pntIndx(:);
    numPoints = length(pntIndx);
    
    numCurves = size(conn,1);
    
    % find line loops and tria mesh used to define surface (Gmsh)
    [ phaseLoops, phaseTria ] = tria2Surface( vert,conn,ele,tnum );
    
    numSurfaces = 0;
    for i = 1: length(phaseLoops)
        numSurfaces = numSurfaces + length(phaseLoops{i});
    end
    
    numVolumes = 0;
    
    % ---------------------------------------------------------------------
    fprintf( fid, '$Entities\n' );
    fprintf( fid, [...
        '%d %d %d %d'  '\n'...
        ], ...
        numPoints, numCurves, numSurfaces, numVolumes );
    
    % ---------------------------------------------------------------------
    % Step2.1. point
    % ---------------------------------------------------------------------
    zeroVec = zeros(numPoints,1);   % column vector
    
    % pnts: NPx4
    pnts = [ pntIndx, vert(pntIndx,:), zeroVec, zeroVec ];

    % pointTag(int) X(double) Y(double) Z(double)
    %   numPhysicalTags(size_t) physicalTag(int)
    % example:
    % 1 0.5 0.5 0 0 
    % 2 0.5 3.5 0 0 
    
    fprintf( fid, ...
        [ ...
        fmtNodeNum, ' ', fmtNodeCo, ' ', fmtNodeCo, ' %d %d', '\n' ...
        ], ...
        pnts' ...
        );

    % end point
    % ---------------------------------------------------------------------
    % Step2.2. curve
    % ---------------------------------------------------------------------
    % Curves are constraining edges (line loop) in my case.
    % conn is the constraining edges

    curveTag = 1: numCurves;    % row
    curveTag = curveTag';       % column vector
    zeroVec = zeros(numCurves,1);   % column vector
    all2Vec = 2 + zeroVec;          % column vector

    % get coordinates of two end point of an boundary edge
    idx1 = conn(:,1);
    idx2 = conn(:,2);
    coor1 = vert(idx1,:);   % NCx2
    coor2 = vert(idx2,:);   % NCx2
    
    % curves: NCx11
    curves = [curveTag, coor1, zeroVec, coor2, zeroVec, ...
                zeroVec, all2Vec, idx1, -idx2 ];

    % curveTag(int) minX(double) minY(double) minZ(double)
    %   maxX(double) maxY(double) maxZ(double)
    %   numPhysicalTags(size_t) physicalTag(int) ...
    %   numBoundingPoints(size_t) pointTag(int; sign encodes orientation)
    % example:
    % 1 0.5 0.5 0 0.5 3.5 0 0 2 1 -2 
    % 2 0.5 0.5 0 2.5 0.5 0 0 2 1 -3 

    fprintf( fid, ...
        [ ...
        '%d ', fmtNodeCo, ' ', fmtNodeCo, ' %d ', fmtNodeCo, ' ', fmtNodeCo, ' %d', ...
        ' %d %d %d %d', '\n' ...
        ], ...
        curves' ...
        );

    % end curve
    % ---------------------------------------------------------------------
    % Step2.3. surface
    % ---------------------------------------------------------------------

    % % surfaces = cell(numSurfaces,1);
    nSurf = 1;  % counter

    % phaseLoops{i} is one phase. phaseLoops{i}{j} is one surface.
    % Build & print surfaces 
    for i = 1: length(phaseLoops)
        for j = 1: length(phaseLoops{i})
            tempSurf = phaseLoops{i}{j};
            outerLineLoopIdx = tempSurf{1};
            
            % extract coordinates
            vertIdx1 = conn( abs(outerLineLoopIdx),1 );
            vertIdx2 = conn( abs(outerLineLoopIdx),2 );
            vert1 = vert(vertIdx1,:);
            vert2 = vert(vertIdx2,:);
            tempVert = [ vert1; vert2 ];
            
            % get extreme
            minXY = min(tempVert);
            maxXY = max(tempVert);
            minX = minXY(1);
            minY = minXY(2);
            maxX = maxXY(1);
            maxY = maxXY(2);
            
            physicalTag = i;    % the i-th phase
            
            curveTag = vertcat(tempSurf{:});    % column vector
            curveTag = curveTag';               % row vector

            % surfaceTag(int) minX(double) minY(double) minZ(double)
            %   maxX(double) maxY(double) maxZ(double)
            %   numPhysicalTags(size_t) physicalTag(int) ...
            %   numBoundingCurves(size_t) curveTag(int; sign encodes orientation)
            % example:
            % 1 2.5 0.5 0 5.5 3.5 0 1 1 4 5 7 -6 -4 
            % 2 0.5 0.5 0 2.5 3.5 0 1 2 4 2 4 -3 -1 
            % 3 5.5 0.5 0 8.5 3.5 0 1 2 4 7 9 -10 -8 

            % surfaces{nSurf} = [ nSurf, minX, minY, 0, ...
            %                     maxX, maxY, 0, ...
            %                     1, physicalTag, length(curveTag), ...
            %                     curveTag ...
            %                     ];

            fprintf( fid, ...
                [ ...
                '%d ', fmtNodeCo, ' ', fmtNodeCo, ' %d ', ...
                fmtNodeCo, ' ', fmtNodeCo, ' %d', ...
                ' %d %d %d ' ...
                ], ...
                nSurf, minX, minY, 0, ...
                maxX, maxY, 0, ...
                1, physicalTag, length(curveTag) ...
                );
            
            fprintf(fid, '%d ', curveTag);
            
            fprintf(fid, '\n');
            
            nSurf = nSurf + 1;
        end
    end

    if nSurf ~= numSurfaces + 1
        error("Somthing wierd happen.")
    end

    % end surface
    % ---------------------------------------------------------------------
    fprintf( fid, '$EndEntities\n' );

    % ---------------------------------------------------------------------
    % Step3. Nodes
    % ---------------------------------------------------------------------
    fprintf( fid, '$Nodes\n' );
    
    % numEntityBlocks(size_t) numNodes(size_t)
    %   minNodeTag(size_t) maxNodeTag(size_t)
    % example:
    % 1 22 1 22

    numNode = size(vert,1);
    fprintf( fid, '%d %d %d %d\n', ...
            1, numNode, 1, numNode );

    % entityDim(int) entityTag(int) parametric(int; 0 or 1)
    %   numNodesInBlock(size_t)
    % example:
    % 2 1 0 22
    
    fprintf( fid, '%d %d %d %d\n', ...
            2, 1, 0, numNode );

    % ---------------------------------------------------------------------
    % nodeTag(size_t)
    % example:
    % 1
    % 2
    % 3

    nodeTag = 1: numNode;
    nodeTag = nodeTag(:);

    fprintf( fid, '%d\n', nodeTag );

    % ---------------------------------------------------------------------
    % node coordinates x y z
    % example:
    % 0.5 0.5 0
    % 0.5 3.5 0

    fprintf( fid, ...
            [ fmtNodeCo, ' ', fmtNodeCo, ' %d', '\n' ], ...
            [ vert, zeros(numNode,1) ]' ...
            );

    fprintf( fid, '$EndNodes\n' );

    % ---------------------------------------------------------------------
    % Step4. Elements
    % ---------------------------------------------------------------------
    fprintf( fid, '$Elements\n' );
    
    % numEntityBlocks(size_t) numElements(size_t)
    %   minElementTag(size_t) maxElementTag(size_t)
    % example:
    % 3 30 1 30

    numEle = size(ele,1);
    fprintf( fid, '%d %d %d %d\n', ...
            numSurfaces, numEle, 1, numEle );


    % entityDim(int) entityTag(int) elementType(int; see below)
    %   numElementsInBlock(size_t)
    %   elementTag(size_t) nodeTag(size_t) ...
    % example:
    % 2 1 2 14
    % 1 12 3 20 
    % 2 2 11 19 
    % 3 13 2 19 

    nSurf = 1;  % surface counter
    nEle = 1;   % element counter
    
    for i = 1: length(phaseTria)
        for j = 1: length(phaseTria{i})
            % surfTria is triangles of one surface (= one entity)
            surfTria = phaseTria{i}{j};
            
            eleType = 2; % 3-node triangle.
            numEleInSurf = size( surfTria, 1 );
            
            fprintf( fid, '%d %d %d %d\n', ...
                    2, nSurf, eleType, numEleInSurf );
            
            % print elements in one surface
            for k = 1: numEleInSurf
                fprintf( fid, '%d %d %d %d\n', ...
                        nEle, surfTria(k,1), surfTria(k,2), surfTria(k,3) );
                nEle = nEle + 1;
            end

            nSurf = nSurf + 1;
        end
    end

    % check
    if nSurf ~= numSurfaces + 1
        error('Bad surface counting');
    end

    if nEle ~= numEle + 1
        error('Bad element counting');
    end

    fprintf( fid, '$EndElements\n' );
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    fclose(fid);
	disp('printMsh Done! Check the msh file!');
    % ---------------------------------------------------------------------

end


