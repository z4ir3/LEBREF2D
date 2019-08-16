function [MESHX] = adjustunstructmesh(MESHX)
%ADJUSTUNSRTUCTMESH renumbering nodes (i.e., edges) for unstructured meshes  
%  output
% --------
%  MESHX : mesh data structure with nodes (edges) renumbered
%          the structure also contains new field for 
%          boundary elements: MESHX.elbnd
%  input
% --------
%  MESHX : mesh data structure 
%
% For a given input unstructured mesh MESHX, this functions changes the order 
% of the nodes per each element in MESH.elem. Basically, it changes the local 
% order of the edges. The reason for this is that the mesh refinement  
% routine assume that the 2nd edge of a given element is the longest
% one and this is not guaranteed by the initial unstructured mesh generated.
%
% NOTE that both positions of nodes and connectivity *do not* change!
%
% ------------------------------------------------------------------
% EXAMPLE: run the example and compare the 2nd, 3rd, 4th row of the element
% mapping matrices before and after the call to the function.
%   MESHX.coord = [0 0.55; 0.4 0; 0.75 0.2; 0.85 0.6; 1 1; 0.35 0.85; 0.5 0.5];
%   MESHX.elem  = [7 6 1; 7 1 2; 7 2 3; 7 3 4; 6 7 4; 4 5 6];
%   MESHX.elem
%   [MESHX] = adjustunstructmesh(MESHX);
%   MESHX.elem
%   plotmesh(MESHX,'',1,1);
% ------------------------------------------------------------------
%
% Function(s) called: detailgrid
%
% LEBREF2D function; Copyright (c) L. Rocchi  
 
  xy  = MESHX.coord;   % elements
  evt = MESHX.elem;    % coordinates
  nel = size(evt,1);   % number of elements
  
% Recover local coordinates
  xlv = zeros(nel,3);
  ylv = zeros(nel,3);
  for ivtx = 1:3
      xlv(:,ivtx) = xy( evt(:,ivtx), 1 ); 
      ylv(:,ivtx) = xy( evt(:,ivtx), 2 );
  end

% -----------------------------------------------------------------------------  
% Compute edge lengths
% -----------------------------------------------------------------------------
  els(:,1) = sqrt( (xlv(:,3) - xlv(:,2)).^2 + (ylv(:,3) - ylv(:,2)).^2 ); % first edge's length
  els(:,2) = sqrt( (xlv(:,1) - xlv(:,3)).^2 + (ylv(:,1) - ylv(:,3)).^2 ); % second edge's length
  els(:,3) = sqrt( (xlv(:,2) - xlv(:,1)).^2 + (ylv(:,2) - ylv(:,1)).^2 ); % third edge's length
  
% -----------------------------------------------------------------------------    
% Find the longest edge per element 
% -----------------------------------------------------------------------------
% This is done by finding the position of maximum length per element:
% - lonedge contains the longest edge per element; 
  [~,lonedge] = max(els,[],2);
  
% -----------------------------------------------------------------------------    
% Changing the order of evt node per element
% -----------------------------------------------------------------------------
  evtrep = [evt, evt(:,1), evt(:,2)];
% - if the longest edge is 1 then take columns [3,4,5] of evtrep;
% - if the longest edge is 2 then take columns [1,2,3] of evtrep, i.e., no changes;
% - if the longest edge is 3 then take columns [2,3,4] of evtrep.
  evt(lonedge==1,:) = evtrep(lonedge==1,[3,4,5]);
  evt(lonedge==2,:) = evtrep(lonedge==2,[1,2,3]);
  evt(lonedge==3,:) = evtrep(lonedge==3,[2,3,4]);
  
% % ---------------------------------------------------------------------------
%   % !DEBUG! Check again the lengths with the new evt: maximum has to be on 2nd position 
%   % Recover local coordinates
%   xlv = zeros(nel,3);  ylv = zeros(nel,3);
%   for ivtx = 1:3
%       xlv(:,ivtx) = xy( evt(:,ivtx), 1 ); 
%       ylv(:,ivtx) = xy( evt(:,ivtx), 2 );
%   end
%   els(:,1) = sqrt( (xlv(:,3) - xlv(:,2)).^2 + (ylv(:,3) - ylv(:,2)).^2 );
%   els(:,2) = sqrt( (xlv(:,1) - xlv(:,3)).^2 + (ylv(:,1) - ylv(:,3)).^2 );
%   els(:,3) = sqrt( (xlv(:,2) - xlv(:,1)).^2 + (ylv(:,2) - ylv(:,1)).^2 );
%   [~,lonedge] = max(els,[],2);
%   if find(lonedge~=2) > 0, error('LEBREF2D: Bad reshaping of evt matrix'); end
% % ---------------------------------------------------------------------------

% Update data structure elements 
  MESHX.elem = evt;
  
% -----------------------------------------------------------------------------    
% Boundary elements 
% -----------------------------------------------------------------------------   
  [MESHY,~]     = detailgrid(MESHX);  
  [belem,bedge] = find( ismember(MESHY.elem ,MESHY.bnd) );
  
% Update data structure boundary elements  
  MESHX.elbnd = [belem, bedge];
  
end % end function  