function [edgenodes] = edge2nodes(MESHX,edgenum,edgelep)
%EDGE2NODES returns the extremal nodes of the edgenum-th of the mesh
%  output
% -----------
%  edgenodes : extremal nodes of the edgenum-th of the mesh
%
%  input 
% ------------
%  MESHX     : mesh data structure
%  edgenum   : edge's number
%  edgelep   : (optional) midpoint element-position matrix
%
% Function(s) called: detailgrid
%
% LEBREF2D function; 12 September 2018
% Copyright (c) 2018 L. Rocchi

  if nargin < 3
      % Generate edgelep
      [~,edgelep] = detailgrid(MESHX);
  end

  edgepos = edgelep(edgenum,3);
  
  if     edgepos == 1, posnodes = [2,3];
  elseif edgepos == 2, posnodes = [3,1];
  elseif edgepos == 3, posnodes = [1,2];
  end
  
  edgenodes = MESHX.elem( edgelep(edgenum,1), posnodes );

end % end function