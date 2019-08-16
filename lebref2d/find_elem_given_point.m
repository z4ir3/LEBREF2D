function [el] = find_elem_given_point(xp,yp,xy,evt)
%FIND_ELEM_GIVEN_POINT returns the element containing a given input point
%
% [el] = find_elem_given_point(xp,yp,xy,evt)
%
% input: 
%       xp,yp    coordinates of the input point
%          xy    vertex coordinate vector 
%         evt    element mapping matrix
%
% output:
%          el    element containg the input point
%     
%    TIFISS function: LR; 22 January 2019
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  nvtx = size(xy,1);    % number of vertices

% Computes distances of the input point (xp,yp) from all vertices of the mesh
  xxp = repmat(xp,nvtx,1);
  yyp = repmat(yp,nvtx,1);
  distances = sqrt( ( xxp - xy(:,1) ).^2 + ( yyp - xy(:,2) ).^2 ); 
  
% Find the vertex of the mesh closer to the input point (xp,yp)  
  [~,vtx] = min( distances );  
% Note that if ~isequal(nvtx,nnz(distances)), 
% then the input point (xp,yp) is a vertex of the mesh, i.e. it is 
% exactly the vtx-th vertex of the mesh
  
% Find the input point-patch, i.e. the elements sharing the vtx-th mesh-point   
  elems = find( sum( (evt==vtx) , 2 ) );

% Check which element of elems the input node (xp,yp) belongs to 
% Using barycentric coordinates for each element in elems:
%
%   xp = a*x1 + b*x2 + c*x3
%   yp = a*y1 + b*y2 + c*y3
%   1  = a + b + c  
%
% where (x1,y1), (x2,y2), and (x3,y3) are the coordinates of the (local) 
% 1st, 2nd, and 3rd vertex of a triangle in elems;
%
% Solution:
%   det := ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
%   a    = ((y2 - y3)*(xp - x3) + (x3 - x2)*(yp - y3)) / det;
%   b    = ((y3 - y1)*(xp - x3) + (x1 - x3)*(yp - y3)) / det;
%   c    = 1 - a - b;
%
% * If 0<=a<=1, 0<=b<=1, and 0<=c<=1 then (xp,yp) belongs to the element.
% Note that this holds also if (xp,yp) belongs to an edge of the element as
% well as if it is a vertex of the element; we take account of 
% this last possibility
  answers = zeros(length(elems),1); 

% Solving for the input point patch
  for elk = 1:length(elems)
      %
      coordel = xy( evt( elems(elk) ,:), : );
      xx = coordel(:,1);
      yy = coordel(:,2);
      %
      % solving the associated linear system
      det = ( yy(2) - yy(3) )*( xx(1) - xx(3) ) + ( xx(3) - xx(2) )*( yy(1) - yy(3) );
      a   = ( (yy(2) - yy(3) )*( xp - xx(3) ) + ( xx(3) - xx(2) )*( yp - yy(3) ) ) / det;
      b   = ( (yy(3) - yy(1) )*( xp - xx(3) ) + ( xx(1) - xx(3) )*( yp - yy(3) ) ) / det;
      c   = 1 - a - b;
      % 
      if ( 0.0 <= a && a <= 1.0 ) && ( 0.0 <= b && b <= 1.0 ) && ( 0.0 <= c && c <= 1.0 )
          % the point (xp,yp) belong to the elems(elk)-th triangle 
          answers(elk) = 1;
      end    
  end

% Find the element which the input node (xp,yp) belongs to
  el = elems( answers==1 );
  
% if length(el) == 1, the input node belong to the el-th element  
  if length(el) > 1
      % = 2, if the input node belongs to the edge shared by the two elements in el
      % = 4, if input node is a vertex of the mesh, i.e. it belongs to all elements in its patch
      el = el(1);
  end

end % end function

% ------------------------------------------------------------------------------
% This is an example of how to use the function to compute the Galerkin
% solution in a point (which is likely not going to be a vertex of the
% mesh) using P2 approximations in the adaptive loop.
% After the solve step, i.e. after that x_gal is computed, the following
% commands would do the job:
%
%       % find the element which the input point belongs to
%       [el] = find_elem_given_point(xp, yp, xy(1:size(xyp1,1),:) , evt(:,1:3));
%       %
%       % get P2 coordinates of the elements      
%       coordelP2 = xy( evt(el,:), : );
%       locx = coordelP2(:,1);
%       locy = coordelP2(:,2);
%       %
%       % Get P2 basis functions on the element
%       BB = [ones(6,1), locx, locy, locx.^2, locx.*locy, locy.^2];
%       psi = zeros(6,1);
%       for k = 1:6
%           b = zeros(6,1);
%           b(k) = 1.0;
%           % coefficients of the psi(k) P2 basis function on the element el
%           a = BB \ b;
%           % psi(k) evaluated at given input point
%           pvalues = [1; xp; yp; xp^2; xp*yp; yp^2];
%           psi(k) = a' * pvalues;
%       end
%       %
%       % P2 Galerkin solution at the given point
%       P2valpoint(refcont) = x_gal( evt(el,:) )' * psi; 
