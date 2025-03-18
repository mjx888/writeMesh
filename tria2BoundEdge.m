function [edge, phaseEdge] = tria2BoundEdge( tria, tnum )
% tria2BoundEdge: convert triangular mesh to boundary edges of surfaces (Gmsh entities).
%
% input:
%   tria: Mesh elements (for linear element). For triangular elements, 
%           it s a Ne-by-3 matrix, where Ne is the number of elements in 
%           the mesh. Each row in tria contains the indices of the nodes 
%           for that mesh element.
%     
%   tnum: Label of phase, which corresponds to physical surface tag in Gmsh. 
%         tnum is a Ne-by-1 array, where Ne is the number of elements.
%         tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%
% output:
%  edge - E-by-2 array. Node numbering of two connecting vertices of
%         boundary edges in all surfaces. Each row is one edge.
%
%
% Copyright (C) 2019-2025 by Jiexian Ma, mjx0799@gmail.com
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % --------------------------------------------------------------------
    % get phaseEdge (nesting cell array)
    num_phase = length( unique(tnum) );
    phaseEdge = cell( 1, num_phase );  % boundary edges of surfaces in phase

    for i = 1: num_phase
	    triaP = tria( tnum==i, : );
	    % triangular mesh of one phase to surface boundary edges
	    phaseEdge{i} = triaPha2BEdge(triaP);
    end

    % --------------------------------------------------------------------
    % convert phaseEdge to array (NE-by-2)
    % Flatten the nested structure into one cell array
    allCells = cellfun(@(x) vertcat(x{:}), phaseEdge, 'UniformOutput', false);
    
    % Concatenate all data into a single NA-by-2 array
    edge = vertcat(allCells{:});

    % remove repeated rows
    edge = unique(edge, 'rows');
    % --------------------------------------------------------------------
end

function edgeCell = triaPha2BEdge(triaP)
% triaPha2BEdge: convert triangular mesh of one phase to surface boundary
% edges
    
    % label isolated mesh regions
    labels = findIsolatedMeshRegions( triaP ) ;
    num_region = length( unique(labels) );
    
    edgeCell = cell( 1, num_region );   % boundary line loop of surface
    
    for i = 1: num_region
        triaIso = triaP( labels == i, : );

        % convert isolate triangular mesh to boundary loops 
        % of a surface region
        edgeCell{i} = findBoundaryEdges(triaIso);
    end
end

function components = findIsolatedMeshRegions( T )
% findIsolatedMeshRegions: Identify connected components where triangles 
% sharing only a vertex are considered disconnected (isolated).
%
% INPUTS:
%   T : (nTri x 3) array of triangle vertex indices
%
% OUTPUT:
%   components : (nTri x 1) array, where components(i) is the index of
%                the connected component to which triangle i belongs,
%                based on sharing *edges* (2 vertices), not just a point.

    nTri = size(T, 1);

    %----------------------------------------------------------------------
    % STEP 1: Build adjacency list so triangles are neighbors only if
    %         they share an entire edge (2 vertices).
    %----------------------------------------------------------------------

    % 1A) Create a map from "edge" -> list of triangle indices that use that edge.
    %     We'll store edges in a canonical form [min, max].
    edgeMap = containers.Map('KeyType','char','ValueType','any');

    for triIdx = 1:nTri
        v = T(triIdx, :);
        % The 3 edges of triangle triIdx:
        edgesTri = [v(1), v(2);
                    v(2), v(3);
                    v(3), v(1)];
        % Sort each row so edge = [min,max]
        edgesTri = sort(edgesTri, 2);

        for e = 1:3
            ePair = edgesTri(e,:);
            key = sprintf('%d-%d', ePair(1), ePair(2));  % e.g. "3-5"
            if ~edgeMap.isKey(key)
                edgeMap(key) = triIdx;  % store just one triangle, or a list
            else
                val = edgeMap(key);
                if ~iscell(val)
                    % convert to cell array if needed
                    val = {val};
                end
                val{end+1} = triIdx; %#ok<AGROW>
                edgeMap(key) = val;
            end
        end
    end

    % 1B) From the edgeMap, build an adjacency list for each triangle.
    adjacency = cell(nTri, 1);
    keysList = edgeMap.keys;
    for i = 1:numel(keysList)
        key = keysList{i};
        val = edgeMap(key);

        if ~iscell(val)
            % Only one triangle had this edge => no adjacency
            continue;
        end

        % If multiple triangles share this same edge, they are neighbors
        triList = [val{:}];  % array of triangle indices
        % Connect all pairs in triList
        for a = 1:length(triList)
            for b = a+1:length(triList)
                tA = triList(a);
                tB = triList(b);
                adjacency{tA}(end+1) = tB; 
                adjacency{tB}(end+1) = tA; 
            end
        end
    end

    % Remove duplicates
    for triIdx = 1:nTri
        adjacency{triIdx} = unique(adjacency{triIdx});
    end

    %----------------------------------------------------------------------
    % STEP 2: Find connected components (based on the adjacency via edges).
    %----------------------------------------------------------------------
    visited = false(nTri, 1);
    components = zeros(nTri, 1);
    currentComp = 0;

    for startTri = 1:nTri
        if ~visited(startTri)
            currentComp = currentComp + 1;
            % Depth-first search (DFS) or BFS
            stack = [startTri];
            visited(startTri) = true;
            components(startTri) = currentComp;

            while ~isempty(stack)
                thisTri = stack(end);
                stack(end) = [];  % pop

                % Neighbors of thisTri
                nbrs = adjacency{thisTri};
                for nb = nbrs
                    if ~visited(nb)
                        visited(nb) = true;
                        components(nb) = currentComp;
                        stack(end+1) = nb; %#ok<AGROW> % push
                    end
                end
            end
        end
    end
end

function boundaryEdges = findBoundaryEdges(Tsub)
% findBoundaryEdges: Identify boundary edges of a sub-mesh.
%
%  INPUT:
%     Tsub  : (nTri x 3) array of triangle vertex indices
%             describing one connected sub-mesh.
%  OUTPUT:
%     boundaryEdges : (nBdryEdges x 2) array of vertex indices
%                     for edges on the boundary.
%                     Each row is [v1, v2] with v1 < v2.
%
%  NOTE:
%  - If your sub-mesh is a subset of a larger mesh, ensure Tsub uses
%    the same vertex indices as the original mesh. Coordinates are
%    not needed in order to find boundary edges.

    nTri = size(Tsub, 1);
    
    % Pre-allocate space for all edges: each triangle has 3 edges
    allEdges = zeros(3*nTri, 2);
    
    % Fill array of edges
    idx = 1;
    for i = 1:nTri
        v = Tsub(i,:);
        % Triangular face => edges: (v1,v2), (v2,v3), (v3,v1)
        triEdges = [v(1),v(2);
                    v(2),v(3);
                    v(3),v(1)];
        
        % Sort each edge so edge = [min, max]
        triEdges = sort(triEdges, 2);
        
        % Store
        allEdges(idx:idx+2, :) = triEdges;
        idx = idx + 3;
    end
    
    % Find unique edges and count occurrences
    [uniqueEdges, ~, ic] = unique(allEdges, 'rows', 'stable');
    counts = accumarray(ic, 1);
    
    % Boundary edges are those that appear exactly once
    isBoundary = (counts == 1);
    boundaryEdges = uniqueEdges(isBoundary, :);
end