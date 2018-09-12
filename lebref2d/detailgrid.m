function [MESHY,edgelep] = detailgrid(MESHX,ixy)
%DETAILGRID linear detail space Y grid generator
%  output 
% ----------
%  MESHY   : mesh data structure for detail grid
%  edgelep : midpoint element-position matrix
%
%  input 
% ----------
%  MESHX   : mesh data structure for the current mesh 
%  ixy     : (optional) flag to recover the xyY coordinates (1/0)
%
% The function creates the data structure MESHY of the detail grid 
% the input mesh, i.e., the structure containg information about 
% all midpoints.
%
% The edgelep matrix contains the elements sharing the 
% midpoints of the mesh as well as the corresponding edge-position.                   
%        
% LEBREF2D function; 12 September 2018
% Copyright (c) 2018 L. Rocchi
  
  if nargin < 2
      % do not run STEP 4 for xyY coordinates; see below
      ixy = 1;
      if nargin < 1
          error('LEBREF2D: at least one input required!');
      end
  end

  elem = MESHX.elem;
  xy   = MESHX.coord;
  nel  = size(elem,1);  % number of elements
  nvtx = size(xy,1);    % number of vertices

% -------------------------------------------------------------------------
% STEP 1: create the edgelep matrix (see description below)
% -------------------------------------------------------------------------
  
% Create the adjacency matrix of pairs of elements sharing the same edge  
% Ex:    (i,j) = 3   and   (j,i) = 5
% It means that elements 3 and 5 share the edge i-j (or j-i). The matrix
% has "symmetric" positions since edge i->j is, for example, the 1st edge
% of element 3, and j->i is the 3rd edge of element 5.
  adj = sparse(nvtx,nvtx);
  adj = adj + sparse(elem(:,2), elem(:,3), 1:nel, nvtx,nvtx);   % all first edges
  adj = adj + sparse(elem(:,3), elem(:,1), 1:nel, nvtx,nvtx);   % all second edges
  adj = adj + sparse(elem(:,1), elem(:,2), 1:nel, nvtx,nvtx);   % all third edges
 
% Create the element-connectivity-matrix eex 
% Ex:    i-th row    ->    14    6    9
% It means that the neighbours of the i-th element are elements 14, 6, and 9, 
% and they lie on the 1st, 2nd, and 3rd edge, respectively (they are "sorted" 
% with respect to edges)
% Ex:    i-th row    ->    14    0    9
% It means that the i-th element is a boundary element and the boundary edge
% is the 2nd one.
  idx_one    = {elem(:,3), elem(:,2)};
  ind_two    = {elem(:,1), elem(:,3)};
  ind_three  = {elem(:,2), elem(:,1)};
  elem_one   = full( adj( sub2ind(size(adj),idx_one{:}) ) );  
  elem_two   = full( adj( sub2ind(size(adj),ind_two{:}) ) );   
  elem_three = full( adj( sub2ind(size(adj),ind_three{:}) ) );
  eex        = [elem_one, elem_two, elem_three];
  
% Replace zeros with the same boundary elements of the corresponding rows:
% this is needed in order to create the next matrix edloc
  [ii,jj] = find(eex == 0);
  idx     = {ii, jj};
  eex( sub2ind(size(eex),idx{:}) ) = ii;
  
% Create the adjacency edge-location matrix edloc (nel-by-nel sparse)
% -------------------------------------------------------------------------
% Ex:    i-th row  ->   0   0   (i,p)=3   0   0   (i,q)=1   0   (i,r)=2
% It means that the neighbours of the i-th element are elements p,q,r and 
% - element p lies on the 3rd edge of i;
% - element q lies on the 1st edge of i;
% - element r lies on the 2nd edge of i.
% -------------------------------------------------------------------------
% If we look at edloc column-wise...
% Ex:    j-th col ->   0   0   (p,j)=1   0   0   (q,j)=3   0   (r,j)=2
% It means that the neighbours of j-th element are elements p,q,r and
% - element j lies on the 1st edge of element p;
% - element j lies on the 3rd edge of element q;
% - element j lies on the 2nd edge of element r.
  edloc = sparse(nel,nel);
  for j = 1:3
      edloc = edloc + sparse((1:nel)', eex(:,j), j*ones(nel,1), nel,nel);
  end
% The boundary elements have the boundary edges along the main diagonal
  
% Take the upper part of edloc excluding the main diagonal (boundary edges)
  uppPart = triu(edloc,1);
  [elem1,elem2,pos1] = find(uppPart);
  
% Find the position with respect to elements in elem2 by taking the entries  
% in the lower symmetric positions of edloc
  idx  = {elem2, elem1};
  pos2 = full( edloc( sub2ind(size(edloc),idx{:}) ) ); 
  
% Final edgelep matrix (#internaledges-by-4 matrix)
% -------------------------------------------------------------------------
% Ex:    i-th row    ->    (i)    7    6    3    2
% It means that the i-th midpoint (Y-basis) will be shared by elements 7 and 6,
% and it will lie on the 3rd edge of element 7 and on the 2nd edge of element 6. 
% Note that at this point edgelep does not include the boundary Y-basis
% functions. It will updated later...
  edgelep = [elem1, elem2, pos1, pos2];
   
% -------------------------------------------------------------------------
% STEP 2: create the uniform refinement mapping matrix elemY
% -------------------------------------------------------------------------
% Ex:    i-th row    ->     8    3    14
% It means that the Y-basis functions (midpoint) of the i-th element are 
% 8, 3, 14, and they lie on the 1st, 2nd, 3rd edge of the element i, respectively.
% It is like MESH.elem for current mesh's vertices, but for midpoints.
% Ex:    i-th row    ->     8    0    14
% It means that the 2nd Y-basis function (midpoint) will be on the boundary
  nYb  = size(edgelep,1); 
  evtY = zeros(nel,3);
  idx1 = {elem1, pos1};
  idx2 = {elem2, pos2};
  evtY( sub2ind(size(evtY), idx1{:}) ) = (1:nYb);
  evtY( sub2ind(size(evtY), idx2{:}) ) = (1:nYb); 
  
% -------------------------------------------------------------------------
% STEP 3: update the edgelep and evtY matrices
% -------------------------------------------------------------------------
% Including the boundary midpoints into edgelep:
% - in evtY, they have to be inserted where currently there are zeros:
%   we insert the new midpoints starting from the last internal midpoint;
% - then, we update edgelep by appending the new rows for the boundary 
%   midpoints: it will be a #alledges-by-4 matrix.

% Number of internal Y-basis 
  nIntBasY = size(edgelep,1);   
  
% Boundary elements
  [bel,beledge] = find(evtY == 0);
  idx = {bel, beledge};

% Update evtY
  boundY = nIntBasY+1 : nIntBasY+length(bel);
  evtY( sub2ind(size(evtY),idx{:}) ) = boundY; 
 
% Create new rows for the edgelep matrix corresponding to the boundary midpoints
  newrows = [bel, bel, beledge, beledge]; 

% Appending rows
  edgelep = [edgelep; newrows];

% -------------------------------------------------------------------------
% STEP 4: create the xyY matrix for the coordinates of the midpoints
% -------------------------------------------------------------------------
  if ixy       
      x = xy(:,1);
      y = xy(:,2);
  
      % Extract the midpoints according to their edge-position (1,2 or 3) 
      % wrt the elements in the first column of edgelep
      [elEdge1,~] = find(edgelep(:,3) == 1);
      [elEdge2,~] = find(edgelep(:,3) == 2);
      [elEdge3,~] = find(edgelep(:,3) == 3);

      % Recovering all the 1st, 2nd, 3rd edges
      firstEd  = elem( edgelep(elEdge1,1), [2,3] );   % first edges
      secondEd = elem( edgelep(elEdge2,1), [3,1] );   % second edges
      thirdEd  = elem( edgelep(elEdge3,1), [1,2] );   % third edges

      % Check that size(frstEd,1) + size(scndEdEd,1) + size(thirdEd,1) == size(YInfo,1)
      if ~isequal( size(firstEd,1) + size(secondEd,1) + size(thirdEd,1), size(edgelep,1))
          error('LEBREF2D:Something wrong with the xyY matrix; check dimensions!');
      end
      
      % Compute midpoint's coordinates
      coordEdge1(:,1) = 0.5 * ( x(firstEd(:,1))  + x(firstEd(:,2))  );    
      coordEdge1(:,2) = 0.5 * ( y(firstEd(:,1))  + y(firstEd(:,2))  );
      %
      coordEdge2(:,1) = 0.5 * ( x(secondEd(:,1)) + x(secondEd(:,2)) );    
      coordEdge2(:,2) = 0.5 * ( y(secondEd(:,1)) + y(secondEd(:,2)) );
      %
      coordEdge3(:,1) = 0.5 * ( x(thirdEd(:,1))  + x(thirdEd(:,2))  );    
      coordEdge3(:,2) = 0.5 * ( y(thirdEd(:,1))  + y(thirdEd(:,2))  );
  
      % Finally create the xy coordinate matrix for the midpoints
      xyY(elEdge1,:) = coordEdge1;
      xyY(elEdge2,:) = coordEdge2;
      xyY(elEdge3,:) = coordEdge3;
      
      % Save xyY coordinates
      MESHY.coord = xyY;
  end
 
% Save structure
  MESHY.elem = evtY;
  MESHY.bnd  = boundY;

end % end function