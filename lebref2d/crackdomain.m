function [MESH] = crackdomain(reflev,slith)
%CRACKDOMAIN structured crack domain (crack on the left) grid generator
% Crack domain generated: (-1,1)^2 \ (-1,0)x{0}.
%  output
% ---------
%  MESH   : mesh data structure
%
%  input 
% ---------
%  reflev : (optional) initial refinement level (default =0)
%  slith  : (optional) half of the height of the slit at (-1,0)x{0}
%
% The slit is created by defining two points:
% - P1 = (xP1,yP1) = (-1,-slith) 
% - P2 = (xP2,yP2) = (-1,+slith)
% which are connected with the "fixed" point Pf=(xf,yf)=(0,0). Then: 
% - smaller values of slitheight -> smaller slits
% - larger values of slitheight  -> large slits 
%
% Default value for slith = 0.01;
%
% -----------------------------------
% EXAMPLES usage:
% - crackdomain;            % default grid and slith
% - crackdomain(1);         % 1 refinement, default slith
% - crackdomain(0,0.03);    % default grid, slith=0.03
% - crackdomain(1,0.03);    % 1 refinement, slith=0.03
% -----------------------------------
%
% See also SQUAREDOMAIN, LSHAPEDOMAIN
%
% LEBREF2D function; Copyright (c) L. Rocchi  
  
  if nargin < 2
      % default value 
      slith = 0.01;
      if nargin < 1
          reflev = 0;
      end
  end

% Fixed point of the slit Pf=(0,0) (the origin)  
  xf = 0.0;
  yf = 0.0;  
  
% ------------------------------------------------------------
% Line connecting P1 with Pf
% ------------------------------------------------------------
% Point P1=(-1,-slith) on the lower left vertical boundary 
  xP1 = -1.0;
  yP1 = -slith;
  
% ------------------------------------------------------------
% Line connecting P2 with Pf
% ------------------------------------------------------------
% Point P2 on the upper left vertical boundary 
  xP2 = -1.0;
  yP2 = +slith;
  
% ------------------------------------------------------------
% Mesh
% ------------------------------------------------------------
% Vertices coordinates map
  xymap = [-1.00    -1.00; ...
            0.00    -1.00; ...
            1.00    -1.00; ...
             xP1      yP1; ...
              xf       yf; ... 
            1.00     0.00; ...
             xP2      yP2; ...
           -1.00     1.00; ...
            0.00     1.00; ...
            1.00     1.00; ...
           -0.50    -0.50; ...
           +0.50    -0.50; ...
           -0.50     0.50; ...
            0.50     0.50]; 
    
% Elements map
  elmap = [ 1    11     4; ...
            2    11     1; ...
            5    11     2; ...
            4    11     5; ...
            2    12     5; ...
            3    12     2; ...
            6    12     3; ...
            5    12     6; ...
            7    13     8; ...
            5    13     7; ...
            9    13     5; ...
            8    13     9; ...
            5    14     9; ...
            6    14     5; ...
           10    14     6; ...
            9    14    10];

% Boundary elements/edges
  elbdmap = [  1    2; ...
               2    2; ...
               4    2; ...
               6    2; ...
               7    2; ...
               9    2; ...
              10    2; ...
              12    2; ...
              15    2; ...
              16    2];
          
% Interior vertices
  intvtx = [11 12 13 14]';
  
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