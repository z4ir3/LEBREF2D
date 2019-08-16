function [MESHX,MESHY,MMele,MMedge,edgelep] = lebmeshref(MESHX,MESHY,edgelep,Mset,imark)  
%LEBMESHREF mesh refinement based on longest edge bisection (LEB) algorithm
%  output 
% ---------
%  MESHX   : mesh data structure (refined mesh)
%  MESHY   : mesh data structure for detail grid of the refined mesh
%  MMele   : set of overall marked elements that are refined
%  MMedge  : set of overall marked edges that are bisected 
%  edgelep : midpoint element-position matrix of the refined mesh
%
%  input 
% ---------
%  MESHX   : mesh data structure (before refinement)
%  MESHY   : mesh data structure for detail grid (before refinement)
%  edgelep : midpoint element-position matrix of the current mesh
%  Mset    : set of marked elements/edges
%  imark   : 1 or 2 if Mset is a vector of elements or edges, respectively
%
% Mesh refinement based on the longest edge bisection (LEB) as a 
% particular case of the newest vertex bisection (NVB) algorithm. 
% The longest edges of marked elements are marked first (reference edges). 
%
% Function(s) called: getallmarkelem
%                     bisection
%                     detailgrid
%                     
% LEBREF2D function; Copyright (c) L. Rocchi  

% Check that Mset is a column vector
  if size(Mset,2) > size(Mset,1)
      Mset = Mset';
  end

% Get the overall set of marked elements/edges in order to keep conformity
  [MMele,MMedge] = getallmarkelem(Mset,MESHY.elem,edgelep,imark); 

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
% Element boundary mapping matrix and new detail grid
% -----------------------------------------------------------------
% Note that, in principle, the mesh-refinement is already finished 
% at this point. Here, we compute the new element boundary mapping 
% matrix and update the detail grid that is used the next iteration
% of the adaptive loop
  [MESHY,edgelep] = detailgrid(MESHX);
  [belem,bedge]   = find( ismember(MESHY.elem, MESHY.bnd) );
  MESHX.elbnd     = [belem, bedge];
 
end  % end function