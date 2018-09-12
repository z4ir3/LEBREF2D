function plotmarkedelem(MESH,Mele,plotelem,plotvtx)
%PLOTMARKEDELEM plots the marked elements of a given mesh 
%  input 
% -----------
%  MESH     : mesh data structure
%  Mele     : vector of marked elements 
%  plotelem : (optional) write elements' number (1/0)
%  plotvtx  : (optional) write vertices' number (1/0)
%
% Default colour is RGB orange [255 140 0]./255;
%
% ----------------------------------------------------------
% EXAMPLE 1. Calling without optional inputs:
%  plotmarkedelem(MESH,Mele);  
%
% EXAMPLE 2:
%  [MESH] = lshapedomain;
%  Mele   = randperm(size(MESH.elem,1),ceil(0.3*size(MESH.elem,1)));
%  plotmarkedelem(MESH,Mele,1,1);  
% ----------------------------------------------------------
%
% See also PLOTMESH, PLOTMESHXANDY
%
% LEBREF2D function; 12 September 2018
% Copyright (c) 2018 L. Rocchi

  if nargin < 4
      % no vertices's number
      plotvtx = '';        
      if nargin < 3
          % no elements' numbers
          plotelem = 0;
          if nargin < 2
              error('LEBREF2D: too less input arguments!');
          end
      end
  end
  
  evt  = MESH.elem;
  xy   = MESH.coord;
  nel  = size(evt,1);   % number of elements
  nvtx = size(xy,1);    % number of vertices
  
% List of available colours
  darkBlue   = [0   141 226]./255;
  darkGreen  = [0.0 100 0.0]./255;
  lightGrenn = [0.0 153 0.0]./255;
  brownPeru  = [205 133  63]./255;
  chocolate  = [210 105  30]./255;
  orange     = [255 140   0]./255;
  markcolor  = orange;
  
% -------------------------------------------------------------  
% Adjacency matrix
% -------------------------------------------------------------  
% This is a sparse adjacency matrix for the mesh given by evt matrix
  adj = sparse(nvtx,nvtx);
  evtt = [evt,evt(:,1)];
  for j = 1:3
      ncol1 = evtt(:,j);
      ncol2 = evtt(:,j+1);
      adj = adj + sparse(ncol1,ncol2,1,nvtx,nvtx);
  end  

% -------------------------------------------------------------    
% Plot mesh
% -------------------------------------------------------------  
  figure;
  gplot(adj,xy,'b');
  hold on;

% Plot coloured marked elements
  LM = length(Mele);
  XY = zeros(3,2,LM);
  for i = 1:LM
      elem      = Mele(i);
      vtxs      = evt(elem,:);
      xyelem    = xy(vtxs,:);
      XY(:,:,i) = xyelem;
      fill(XY(:,1,i),XY(:,2,i),markcolor);
  end
  
% Display elements' numbers  
  if plotelem
      xl_v = zeros(nel,3); 
      yl_v = zeros(nel,3); 
      for ivtx = 1:3
          xl_v(:,ivtx) = xy(evt(:,ivtx),1); % x-coordinates of all vertices
          yl_v(:,ivtx) = xy(evt(:,ivtx),2); % y-coordinates of all vertices
      end
      % Element's centroid coordinates
      xyc(:,1) = (1/3) * sum(xl_v,2);
      xyc(:,2) = (1/3) * sum(yl_v,2);
      %
      % Write element's number
      elem   = (1:nel)';
      elenum = int2str(elem);      
      text(xyc(:,1),xyc(:,2),elenum,'Color','blue','Fontsize',12);
  end
   
% Display vertices' numbers
  if plotvtx      
      vertices = (1:nvtx)';
      vtxnum   = int2str(vertices);
      text(xy(:,1),xy(:,2),vtxnum,'Color','black','Fontsize',14);
  end 
  
  axis('square');
  set(gca,'FontSize',17);
  set(gcf,'units','normalized','Position',[0.25 0.05 0.5 0.7]);
  axis off;
  hold off;
   
end % end function