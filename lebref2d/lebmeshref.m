function [MESHX,MMele] = lebmeshref(MESHX,Mele)  
%LEBMESHREF mesh refinement based on longest edge bisection (LEB) algorithm
%  output 
% ---------
%  MESH   : mesh data structure (refined mesh)
%  MMele  : overall set of marked elements that are refined 
%
%  input 
% ---------
%  MESHX  : mesh data structure (before refinement)
%  Mele   : set of marked elements
%
% Mesh refinement based on the longest edge bisection (LEB) as a 
% particulat case of the newest vertex bisection (NVB) algorithm. 
% The longest edges of marked elements are marked first (reference edges). 
%
% Function(s) called: detailgrid
%                     getallmarkelem
%                     bisection
%                     
% LEBREF2D function; 12 September 2018
% Copyright (c) 2018 L. Rocchi  

% Check that Mele is a column vector
  if size(Mele,2) > size(Mele,1)
      Mele = Mele';
  end

% Detail grid associated with the current mesh
  [MESHY,edgelep] = detailgrid(MESHX);

% Get the overall set of marked elements/edges in order to keep conformity
  [MMele,MMedge] = getallmarkelem(Mele,MESHY.elem,edgelep);                %(Mele,edgelep,MESHY.elem); 

% Overall number of edges
  nedg = size(MESHY.coord,1);
    
% Global numbers of the new nodes (midpoints) to be inserted in the mesh
  newnodes = size(MESHX.coord,1) + (1:length(MMedge));
  
% Midpoints' global numbers: the following vector markedge is nedg-by-1
% sparse and contains either 0 in the i-th position if the i-th edge 
% has not been marked or the (global) number of the midpoint that would 
% be inserted on that edge
  markedge = sparse(MMedge,1,newnodes,nedg,1);
    
% -----------------------------------------------------------------
% Refine the marked elements/edges and create the new element matrix  
% -----------------------------------------------------------------
  [refelem] = bisection(MMele,MMedge,markedge,MESHX.elem,MESHY.elem);
  
% -----------------------------------------------------------------
% New coordinates vector: appending new midpoints' coordinates
% -----------------------------------------------------------------
  refxy = [MESHX.coord; MESHY.coord(MMedge,:)];
  
% -----------------------------------------------------------------
% New interior/boundary nodes
% -----------------------------------------------------------------
% New boundary nodes
  refbnd = [MESHX.bnd; nonzeros( markedge(MESHY.bnd) )];
  
% New interior nodes vector
  totvtx = 1:size(MESHX.coord,1);
  refint = totvtx( ~ismember(totvtx,refbnd) )';
  
% Save the structure
  MESHX.coord = refxy;
  MESHX.elem  = refelem;
  MESHX.int   = refint;
  MESHX.bnd   = refbnd;
  
% -----------------------------------------------------------------
% Element boundary mapping matrix for the refined mesh
% -----------------------------------------------------------------
% This is indeed not needed for the real mesh refinement which is 
% already finished at this stage; 
% It is needed if the new element boundary map is required after 
% mesh refinement.
% In the latter case, the function DETAILGRID is used although it 
% does not need to run up to its final STEP 4 (see inside) which recovers 
% the xyY-coordinates matrix; calling is done with optional input '0';
  [MESHY,~]     = detailgrid(MESHX,0);
  [belem,bedge] = find( ismember(MESHY.elem, MESHY.bnd) );
  MESHX.elbnd   = [belem, bedge];
 
end  % end function