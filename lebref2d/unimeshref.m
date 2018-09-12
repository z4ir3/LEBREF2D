function [MESHREF] = unimeshref(MESH,reftype,iplot)
%UNIMESHREF uniformly refines the triangular mesh using either red or bisec3 refinement
%  output
% ----------
%  MESHREF : mesh data structure (after uniform refinement)
%
%  input
% ----------
%  MESH    : mesh data structure (before uniform refinement)
%  reftype : red or bisec3 (1/2)
%  iplot   : (optional) grid plotting switch
%
% ---------------------------------------------------------
% EXAMPLE 1: red refinement 
%  [MESH] = lshapedomain; 
%  plotmesh(MESH,'',1,1);
%  [MESHREF] = unimeshref(MESH,1);
%  plotmesh(MESHREF,'',1,1);
%
% EXAMPLE 2: bisec3 refinement 
%  [MESH] = lshapedomain; 
%  plotmesh(MESH,'',1,1);
%  [MESHREF] = unimeshref(MESH,2);
%  plotmesh(MESHREF,'',1,1);
% ---------------------------------------------------------
%
% LEBREF2D function; 12 September 2018
% Copyright (c) 2018 L. Rocchi

  if nargin < 3
      iplot = 0;
      if nargin < 2
          % default bisec3 uniform refinement
          reftype = 2;
      end
  end

  nvtx  = size(MESH.coord,1);  % number of vertices of MESH      
  nel   = size(MESH.elem,1);   % number of elements of MESH      
  
% -----------------------------------------------------------------
% STEP 1: getting edges
% -----------------------------------------------------------------

% Total edges (with repetitions)  
  edges = [MESH.elem(:,[2,3]); ...
           MESH.elem(:,[3,1]); ... 
           MESH.elem(:,[1,2])];
  
% Total edges of the triangulation without repetitions
  [edges,~,ic] = unique(sort(edges,2),'rows');    

% Edges per element, i.e., the new nodes (midpoints) per element
  edgenum = [ic(1:nel), ic(nel+1:2*nel), ic(2*nel+1:3*nel)];
    
% -----------------------------------------------------------------
% STEP 2: elements' nodes and midpoints map
% -----------------------------------------------------------------
% This is a nel-by-6 matrix containing MESH.elem on the first 3 columns
% and the associated midpoints on the last 3 columns
  p2nodes = edgenum + nvtx;
  p2evt   = [MESH.elem, p2nodes];
  
% -----------------------------------------------------------------
% STEP 3: nodes' coordinates of the refined mesh
% -----------------------------------------------------------------  
% This is a #allnodes-by-2 matrix constisting of MESH.coord and 
% the (new) midpoints' coordinates appended at the end 
  mxy(:,1) = 0.5 * (MESH.coord(edges(:,1),1) + MESH.coord(edges(:,2),1));
  mxy(:,2) = 0.5 * (MESH.coord(edges(:,1),2) + MESH.coord(edges(:,2),2));
  MESHREF.coord = [MESH.coord; mxy];
  
% -----------------------------------------------------------------
% STEP 4: boundary nodes of the refined mesh
% -----------------------------------------------------------------
% This is a #allboundarynodes vector constisting of MESH.bnd and 
% the boundary midpoints appended at the end
  checkbd     = sum( ismember(edges,MESH.bnd) , 2);
  p2nodesbd   = find(checkbd==2) + nvtx;  
  MESHREF.bnd = [MESH.bnd; p2nodesbd];

% -----------------------------------------------------------------
% STEP 5: element map and boundary element map
% -----------------------------------------------------------------
  
% Allocate memory 
  revt    = zeros(4*nel,3);
  eboundt = MESH.elbnd;
  
  if reftype == 1
      % -----------------------------------------------------------
      % Uniform RED refinement
      % -----------------------------------------------------------
      revt(1:4:(4*nel),:) = p2evt(:,[1,6,5]);
      revt(2:4:(4*nel),:) = p2evt(:,[6,2,4]);  
      revt(3:4:(4*nel),:) = p2evt(:,[5,4,3]);
      revt(4:4:(4*nel),:) = p2evt(:,[4,5,6]);
      %
      % Data structures for the boundary element map of the refined mesh
      %
      % Extract boundary elements according to the boundary edge's number
      bedgeOne   = eboundt(eboundt(:,2) == 1,:);
      bedgeTwo   = eboundt(eboundt(:,2) == 2,:);
      bedgeThree = eboundt(eboundt(:,2) == 3,:);
      %
      % For a given (old) boundary element in the eboundt matrix, 
      % there will be only 2 (out of the 4) new boundary element-children: 
      % according to the old element's edge number, 1, 2, or 3, the 
      % boundary children will be the [2nd, 3rd], [1st, 3rd] or [1st 2nd] 
      % elements created, respectively
      matElOne   = [4*(bedgeOne(:,1)-1),   4*(bedgeOne(:,1)-1)]   + repmat([2 3],size(bedgeOne,1),1);
      matElTwo   = [4*(bedgeTwo(:,1)-1),   4*(bedgeTwo(:,1)-1)]   + repmat([1 3],size(bedgeTwo,1),1);
      matElThree = [4*(bedgeThree(:,1)-1), 4*(bedgeThree(:,1)-1)] + repmat([1 2],size(bedgeThree,1),1);
      %
      % Adding the corresponding edge number per each child boundary
      % element: according to the old element's edge number, 1, 2, or 3, 
      % the boundary edges will be [1 1], [2 2], or [3 3] 
      matElOne   = [matElOne(:,1),   ones(size(matElOne,1),1),       matElOne(:,2),   ones(size(matElOne,1),1)];
      matElTwo   = [matElTwo(:,1),   repmat(2,size(matElTwo,1),1),   matElTwo(:,2),   repmat(2,size(matElTwo,1),1)];
      matElThree = [matElThree(:,1), repmat(3,size(matElThree,1),1), matElThree(:,2), repmat(3,size(matElThree,1),1)];
      
  else% refine_type == 2     
      % -----------------------------------------------------------
      % Uniform BISEC3 refinement
      % -----------------------------------------------------------
      revt(1:4:(4*nel),:) = p2evt(:,[1,6,5]);
      revt(2:4:(4*nel),:) = p2evt(:,[5,6,2]);
      revt(3:4:(4*nel),:) = p2evt(:,[2,4,5]);
      revt(4:4:(4*nel),:) = p2evt(:,[5,4,3]);
      %
      % Data structures for the boundary element map of the refined mesh
      %
      % Extract boundary elements according to the boundary edge's number
      bedgeOne   = eboundt(eboundt(:,2) == 1,:);
      bedgeTwo   = eboundt(eboundt(:,2) == 2,:);
      bedgeThree = eboundt(eboundt(:,2) == 3,:);
      %
      % For a given (old) boundary element in the eboundt matrix, 
      % there will be only 2 (out of the 4) new boundary element-children: 
      % according to the old element's edge number, 1, 2, or 3, the 
      % boundary children will be the [3rd 4th], [1st, 4th] or [1st 2nd] 
      % elements created, respectively
      matElOne   = [4*(bedgeOne(:,1)-1),   4*(bedgeOne(:,1)-1)]   + repmat([3 4],size(bedgeOne,1),1);
      matElTwo   = [4*(bedgeTwo(:,1)-1),   4*(bedgeTwo(:,1)-1)]   + repmat([1 4],size(bedgeTwo,1),1);
      matElThree = [4*(bedgeThree(:,1)-1), 4*(bedgeThree(:,1)-1)] + repmat([1 2],size(bedgeThree,1),1);
      %
      % Adding the corresponding edge number per each child boundary
      % element: according to the old element's edge number, 1, 2, or 3, 
      % the boundary edges will be [3 1], [2 2], or [3 1] 
      matElOne   = [matElOne(:,1),   repmat(3,size(matElOne,1),1),   matElOne(:,2),   ones(size(matElOne,1),1)];
      matElTwo   = [matElTwo(:,1),   repmat(2,size(matElTwo,1),1),   matElTwo(:,2),   repmat(2,size(matElTwo,1),1)];
      matElThree = [matElThree(:,1), repmat(3,size(matElThree,1),1), matElThree(:,2), ones(size(matElThree,1),1)];         
  end

% Elements of the refined mesh  
  MESHREF.elem  = revt;
 
% Boundary element map for the refined mesh
  MESHREF.elbnd = [ matElOne(:,[1,2]);   matElOne(:,[3,4]); ...
                    matElTwo(:,[1,2]);   matElTwo(:,[3,4]); ...
                    matElThree(:,[1,2]); matElThree(:,[3,4] )]; 

% Interior nodes of the refined mesh
  totvtx      = 1:size(MESHREF.coord,1);
  MESHREF.int = totvtx(~ismember(totvtx,MESHREF.bnd));

  if iplot
      % Plot the uniform refined mesh
      plotmesh(MESHREF, [], 1,1);
  end

end % end function