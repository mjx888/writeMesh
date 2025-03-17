function [ phaseLoops, phaseTria ] = tria2Surface( vert,conn,tria,tnum )
% tria2Surface: convert triangular mesh to surfaces (Gmsh entities).
% 
% A surface is defined by line loops or curve loops of boundary edges
%
% input:
%   vert: Mesh nodes (for linear element). It’s a Nn-by-2 matrix, where 
%           Nn is the number of nodes in the mesh. Each row of vert 
%           contains the x, y coordinates for that mesh node.
%     
%   tria: Mesh elements (for linear element). For triangular elements, 
%           it s a Ne-by-3 matrix, where Ne is the number of elements in 
%           the mesh. Each row in tria contains the indices of the nodes 
%           for that mesh element.
%     
%   tnum: Label of phase, which corresponds to physical surface tag in Gmsh. 
%         tnum is a Ne-by-1 array, where Ne is the number of elements.
%         tnum(j,1) = k; means the j-th element belongs to the k-th phase.
%     
%   conn: C-by-2 array of constraining edges, where each row defines an edge
%
% output:
%  phaseLoops = C
%  C is a nesting cell array for storing multiple loops. C is a 1-by-P cell 
%    array. C{i} means the i-th physical surface. C{i} is a 1-by-S cell 
%    array. C{i}{j} means the j-th plane surface within the i-th physical 
%    surface. C{i}{j} is 1-by-L cell array. C{i}{j}{k} means the k-th loop 
%    of the j-th plane surface within the i-th physical surface. C{i}{j}{k} 
%    stores the line indices within a loop. C{i}{j}{k} is an N-by-1 array.
%    line indices are storing in conn.
%
%  phaseTria is a nesting cell array for storing triangular mesh for each
%    surface. phaseTria{i}{j} means the j-th plane surface within the 
%    i-th physical surface. phaseTria{i}{j} is a p-by-3 array.
%
%
% Copyright (C) 2019-2025 by Jiexian Ma
% 
% Project website: https://github.com/mjx888/im2mesh
%

    % convert to loops
    num_phase = length( unique(tnum) );
    phaseLoops = cell( 1, num_phase );  % boundary line loop of surfaces in phase
    phaseTria = cell( 1, num_phase );   % tria mesh of surfaces in phase

    for i = 1: num_phase
	    triaP = tria( tnum==i, : );
	    % triangular mesh of one phase to surface loop cell
	    [ phaseLoops{i}, phaseTria{i} ] = triaPha2loop(triaP, vert, conn);
    end

end

function [surfaceLoops, surfaceTria] = triaPha2loop(triaP, vert, edge)
% triaPha2loop: convert triangular mesh of one phase to surface loops cell
%
    
    % label isolated mesh regions
    labels = findIsolatedMeshRegions( triaP ) ;
    num_region = length( unique(labels) );
    
    surfaceLoops = cell( 1, num_region );   % boundary line loop of surface
    surfaceTria = cell( 1, num_region );    % tria mesh of surface
    
    for i = 1: num_region
        triaIso = triaP( labels == i, : );

        % convert isolate triangular mesh to boundary loops 
        % of a surface region
        surfaceLoops{i} = triaIso2loop( triaIso, vert, edge );
        surfaceTria{i} = triaIso;
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


function loopsEdgesInd = triaIso2loop( triaIso, vert, edge )
% triaIso2loop: convert isolate triangular mesh to boundary loops of a 
%               surface region
% plot loopsEdgesInd by function plotLoopsEdgesInd(loopsEdgesInd, edge, vert);

    boundaryEdges = findBoundaryEdges(triaIso);
    loops = groupBoundaryEdgesIntoLoops(boundaryEdges);
    loops = makeOuterBoundaryFirst(loops, vert);
    % plot loops using function plotLoops(loops, vert);
    
    loopsEdges = convertLoopsToEdgePairs(loops);
    loopsEdgesInd = createLoopsEdgesInd(loopsEdges, edge);
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


function loops = groupBoundaryEdgesIntoLoops(boundaryEdges)
% groupBoundaryEdgesIntoLoops: Group boundary edges into separate loops.
%  
%  boundaryEdges : (nBdryEdges x 2), each row [v1, v2] (no strict requirement of v1 < v2).
%  loops         : cell array of loops; loops{i} is a list of vertex IDs [v1, v2, ..., vK]
%                  that form a closed boundary chain. The last vertex typically equals the first.
%
%  This function attempts to trace out boundary loops edge by edge. It avoids infinite
%  loops by picking one unvisited "next" edge at each vertex. If there's branching
%  (multiple possible next edges), it picks the first one. Any remaining unvisited edges
%  from that branching point will start a new loop in the next iteration.

    if isempty(boundaryEdges)
        loops = {};
        return;
    end

    %--- 1) Build adjacency from vertex -> boundary neighbors
    allVerts = unique(boundaryEdges(:));
    adjMap = containers.Map('KeyType','double','ValueType','any');
    for v = allVerts(:)'
        adjMap(v) = [];
    end
    
    for i = 1:size(boundaryEdges, 1)
        v1 = boundaryEdges(i,1);
        v2 = boundaryEdges(i,2);
        adjMap(v1) = [adjMap(v1), v2];
        adjMap(v2) = [adjMap(v2), v1];
    end
    
    %--- 2) Keep track of visited edges
    nEdges = size(boundaryEdges,1);
    visitedEdges = false(nEdges, 1);

    % We create a map from (v1-v2 or v2-v1) -> edge index,
    % so we can easily find which edge index corresponds to a pair of vertices.
    edgeMap = containers.Map('KeyType','char','ValueType','double');
    for e = 1:nEdges
        v1 = boundaryEdges(e,1);
        v2 = boundaryEdges(e,2);
        k1 = sprintf('%d-%d', v1, v2);
        k2 = sprintf('%d-%d', v2, v1);
        edgeMap(k1) = e;
        edgeMap(k2) = e;  % same edge index for reversed pair
    end
    
    getKey = @(a,b) sprintf('%d-%d', a, b);
    
    %--- 3) Find loops by traversing unvisited edges
    loops = {};
    for e = 1:nEdges
        if ~visitedEdges(e)
            % Start a new loop from this unvisited edge
            v1 = boundaryEdges(e,1);
            v2 = boundaryEdges(e,2);
            
            loop = [v1, v2];
            visitedEdges(e) = true;  % mark this edge visited
            
            currentVertex = v2;
            prevVertex    = v1;
            
            % Keep moving forward around the loop
            while true
                neighbors = adjMap(currentVertex);
                
                % Remove the vertex we came from
                neighbors(neighbors == prevVertex) = [];
                
                % Among these neighbors, find any unvisited edges
                unvisitedNext = [];
                for cand = neighbors
                    eIdx = edgeMap(getKey(currentVertex, cand));
                    if ~visitedEdges(eIdx)
                        unvisitedNext(end+1) = cand; %#ok<AGROW>
                    end
                end
                
                if isempty(unvisitedNext)
                    % No unvisited edge extends from currentVertex
                    % => we can't proceed further to close a loop.
                    break;
                end
                
                % In case of branching, pick the first unvisited neighbor
                nextVertex = unvisitedNext(1);
                
                % Mark that edge visited
                nextEdgeIdx = edgeMap(getKey(currentVertex, nextVertex));
                visitedEdges(nextEdgeIdx) = true;
                
                % Add nextVertex to the loop
                loop(end+1) = nextVertex; %#ok<AGROW>
                
                % Advance
                prevVertex = currentVertex;
                currentVertex = nextVertex;
                
                % Check if we've come full circle
                if currentVertex == v1
                    % Loop is closed
                    break;
                end
            end
            
            loops{end+1} = loop; %#ok<AGROW>
        end
    end
end

function loops = makeOuterBoundaryFirst(loops, P)
% makeOuterBoundaryFirst: Reorder loops so the one with the largest |area| 
% is first.
%
%  loops : cell array, each element is a list of vertex indices defining a boundary loop.
%  P     : (nPts x 2) vertex coordinates.
%
%  This assumes that in your mesh, the loop with the largest absolute area
%  is indeed the outer boundary.

    % 1) Compute signed area of each loop
    [loopArea, ~] = classifyBoundaryLoops(loops, P);
    
    % 2) Find the loop with the maximum absolute area
    [~, mainIdx] = max(abs(loopArea));
    
    % 3) Reorder loops so that this loop is first
    nLoops = numel(loops);
    newOrder = [ mainIdx, setdiff(1:nLoops, mainIdx) ];
    
    loops = loops(newOrder);
end

function [loopArea, loopType] = classifyBoundaryLoops(loops, P)
% classifyBoundaryLoops: Compute signed area of each loop and label as 
% 'outer' or 'hole'.
%
%  loopArea : (nLoops x 1) signed area of each loop
%             > 0 => CCW => "outer"
%             < 0 => CW  => "hole"
%  loopType : cell array, 'outer' or 'hole'

    nLoops = numel(loops);
    loopArea = zeros(nLoops,1);
    loopType = cell(nLoops,1);

    for i = 1:nLoops
        coords = P(loops{i}, :);
        loopArea(i) = polygonSignedArea2D(coords);
        if loopArea(i) > 0
            loopType{i} = 'outer';
        else
            loopType{i} = 'hole';
        end
    end
end

function signedArea = polygonSignedArea2D(coords)
% polygonSignedArea2D: Compute the signed area of a 2D polygon (shoelace 
% formula).
%  Positive => CCW orientation, Negative => CW orientation.

    % Ensure the polygon is closed
    if any(coords(end,:) ~= coords(1,:))
        coords(end+1,:) = coords(1,:);
    end
    
    x = coords(:,1);
    y = coords(:,2);
    signedArea = 0.5 * sum( x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1) );
end


function loopsEdges = convertLoopsToEdgePairs(loops)
% convertLoopsToEdgePairs: Convert each CLOSED loop of vertex indices into 
% Mx2 edge pairs.
%
%  loops : cell array, where each loops{i} = [v1, v2, ..., vN, v1].
%          That is, the loop is already closed (the last vertex repeats the first).
%
%  loopsEdges : cell array of the same size as 'loops'.
%               loopsEdges{i} is an M-by-2 array of vertex indices (each row is one edge).
%               For a loop with length N, you get N-1 edges (the last pair is (vN, v1)).

    nLoops = numel(loops);
    loopsEdges = cell(size(loops));

    for i = 1:nLoops
        thisLoop = loops{i};
        
        % If the loop is closed (v1 == v(end)), then the number of actual edges
        % is length(thisLoop) - 1. Each edge is (v_i, v_{i+1}).
        % For example, [v1, v2, ..., vN, v1] => edges:
        %   (v1, v2)
        %   (v2, v3)
        %   ...
        %   (v_{N-1}, vN)
        %   (vN, v1)
        
        edges = [thisLoop(1:end-1)', thisLoop(2:end)'];
        
        loopsEdges{i} = edges;
    end
end


function loopsEdgesInd = createLoopsEdgesInd(loopsEdges, edges)
% createLoopsEdgesInd: Create a cell array storing the signed edge indices 
% for each loop.
%
%  loopsEdges : cell array where loopsEdges{i} is (M_i x 2), the edges of loop i.
%               Each row is [v1, v2] in the orientation of that loop's boundary.
%
%  edges      : (N x 2) global array of edges.
%               Each row edges(k,:) = [v1, v2] indicates the orientation as well.
%
%  loopsEdgesInd : cell array of the same size as loopsEdges.
%                  loopsEdgesInd{i} is a column vector of length M_i.
%                  For an edge [v1, v2] in loopsEdges{i}, if edges(k,:) = [v1, v2],
%                  then loopsEdgesInd{i}(j) = +k.
%                  If edges(k,:) = [v2, v1], then loopsEdgesInd{i}(j) = -k.
%
%  ASSUMPTION: 1-based indexing for vertices and edges (standard in MATLAB).
%              If you have 0-based elsewhere, adjust accordingly.

    nLoops = numel(loopsEdges);
    loopsEdgesInd = cell(size(loopsEdges));
    
    %---------------------------------
    % 1) Build a map from "v1-v2" to signed edge index
    %---------------------------------
    edgeMap = containers.Map('KeyType','char','ValueType','int32');
    
    nEdgesGlobal = size(edges, 1);
    for k = 1:nEdgesGlobal
        v1 = edges(k,1);
        v2 = edges(k,2);
        
        % orientation as stored
        directKey = sprintf('%d-%d', v1, v2);  
        % reversed orientation
        revKey    = sprintf('%d-%d', v2, v1);
        
        % If someone tries to insert a duplicate key, you'd get an error.
        % For a well-defined mesh, each edge should appear exactly once in a given orientation.
        % But to be safe, you might want to check edgeMap.isKey(directKey) etc.
        
        edgeMap(directKey) = +k;  % same orientation => +k
        edgeMap(revKey)    = -k;  % reversed => -k
    end
    
    %---------------------------------
    % 2) For each loop, build a signed index array
    %---------------------------------
    for i = 1:nLoops
        thisEdges = loopsEdges{i};  % M_i x 2
        M = size(thisEdges,1);
        edgeIdxList = zeros(M,1,'int32');  %# for storing signed indices
        
        for j = 1:M
            v1 = thisEdges(j,1);
            v2 = thisEdges(j,2);
            
            key = sprintf('%d-%d', v1, v2);
            if edgeMap.isKey(key)
                edgeIdxList(j) = edgeMap(key);
            else
                % If we don't find it in the map, then something is off
                % (edge does not exist in the global 'edges' list).
                % You could throw an error or store 0 as a sentinel.
                error('Edge [%d %d] not found in global edges.', v1, v2);
            end
        end
        
        loopsEdgesInd{i} = edgeIdxList;  % store as a column vector
    end
end
