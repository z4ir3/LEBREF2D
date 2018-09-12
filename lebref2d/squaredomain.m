function [MESH] = squaredomain(reflev)
%SQUAREDOMAIN structured square domain grid generator
% Square domain generated: (-1,1)^2
%  output
% ---------
%  MESH   : mesh data structure
%
%  input:
% ---------
%  reflev : (optional) initial refinement level (default =0) 
%
% See also LSHAPEDOMAIN, CRACKDOMAIN
%
% LEBREF2D function; 12 September 2018
% Copyright (c) 2018 L. Rocchi

  if nargin < 1
      reflev = 0;
  end

% ------------------------------------------------------------
% Mesh
% ------------------------------------------------------------
% Vertices coordinates map
  xymap = [-1.00   -1.00; ...
            0.00   -1.00; ...
            1.00   -1.00; ...
           -1.00    0.00; ...
            0.00    0.00; ...
            1.00    0.00; ...
           -1.00    1.00; ...
            0.00    1.00; ...
            1.00    1.00]; 
    
% Elements map
  elmap = [ 5     4     1; ...
            1     2     5; ...
            5     2     3; ...
            3     6     5; ...
            5     6     9; ...
            9     8     5; ...
            5     8     7; ...
            7     4     5];  
 
% Boundary elements/edges
  elbdmap = [  1    1; ...
               2    3; ...
               3    1; ...
               4    3; ...
               5    1; ...
               6    3; ...
               7    1; ...
               8    3];
          
% Interior vertices
  intvtx = [5]';
  
% Boundary vertices
  totvtx = 1:size(xymap,1);
  bdvtx  = totvtx(~ismember(totvtx,intvtx))';
  
% MESH structure
  MESH.coord = xymap;
  MESH.elem  = elmap;
  MESH.int   = intvtx;
  MESH.bnd   = bdvtx;
  MESH.elbnd = elbdmap;  
  
% Uniform refinements (if any)  
  for i = 1:reflev
      [MESH] = unimeshref(MESH,2);
  end
  
end % end function