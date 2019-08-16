function plotmeshxandy(MESHX,MESHY,reftype,plotelem,plotvtx)
%PLOTMESHXANDY plots the mesh and its uniform (red/bisec3) refinement 
%  input 
% -----------
%  MESHX    : mesh data structure for the current mesh
%  MESHY    : mesh data structure for detail grid (midpoints)
%  reftype  : red or bisec3 (1/2) refinement
%  plotelem : (optional) write elements' number (1/0)
%  plotvtx  : (optional) write vertices' number (1/0)
%
% The function plots the current mesh together with the associated 
% midpoint-grid (detail Y space); see also DETAILGRID.
%
% ---------------------------------------------------
% EXAMPLE 1: red subdivision
%  [MESHX]   = squaredomain;
%  [MESHY,~] = detailgrid(MESHX);
%  plotmeshxandy(MESHX,MESHY,1,1,1); 
%
% EXAMPLE 2: bisec3 subdivision
%  [MESHX]   = squaredomain;
%  [MESHY,~] = detailgrid(MESHX);
%  plotmeshxandy(MESHX,MESHY,2,1,1);
% ---------------------------------------------------
%
% See also PLOTMESH, PLOTMARKEDELEM
%
% LEBREF2D function; Copyright (c) L. Rocchi  
 
  if nargin < 5
      % no vertices' numbers
      plotvtx = 0;        
      if nargin < 4
          % no elements' numbers
          plotelem = 0;
          if nargin < 3
              % default bisec3 subdivision
              reftype = 1;
              if nargin < 2
                  error('LEBREF2D: too less input arguments!');
              end  
          end
      end
  end
    
  evt  = MESHX.elem;
  xy   = MESHX.coord;
  evtY = MESHY.elem;
  xyY  = MESHY.coord;
      
  nel   = size(evt,1);      % number of elements
  nvtx  = size(xy,1);       % total number of vertices
  nvtxY = size(xyY,1);      % total number of midpoints (edges)
  
% List of available colours
  redcol     = [255   0   0]./255; % red
  bisec3col  = [255 140   0]./255; % orange   
  darkBlue   = [0   141 226]./255;
  darkGreen  = [0.0 100 0.0]./255;
  lightGrenn = [0.0 153 0.0]./255;
  brownPeru  = [205 133  63]./255;
  chocolate  = [210 105  30]./255;
    
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
  hold on;
  gplot(adj,xy,'b');

% Display element's numbers  
  if plotelem
      xl_v = zeros(nel,3); 
      yl_v = zeros(nel,3); 
      for ivtx = 1:3
          xl_v(:,ivtx) = xy(evt(:,ivtx),1); % x-coordinates of all nodes
          yl_v(:,ivtx) = xy(evt(:,ivtx),2); % y-coordinates of all nodes
      end
      % Element's centroid coordinates
      xyc(:,1) = sum(xl_v,2) / 3;
      xyc(:,2) = sum(yl_v,2) / 3;
      %
      % Write element's number
      elem = (1:nel)';
      elenum = int2str(elem);      
      text(xyc(:,1),xyc(:,2),elenum,'Color','blue','Fontsize',12);
  end

% Display vertices' numbers
  if plotvtx      
      vertices = (1:nvtx)';
      vtxnum = int2str(vertices);
      text(xy(:,1),xy(:,2),vtxnum,'Color','black','Fontsize',14);
  end   

% -------------------------------------------------------------------
% Plot Y space w.r.t. the uniform subdivision
% -------------------------------------------------------------------   
  vtxY   = (1:nvtxY)';
  vtxnum = int2str(vtxY);
  if reftype == 1
      %
      % Red subdivision
      %
      % Write vertices' numbers
      text(xyY(:,1),xyY(:,2),vtxnum,'Color',redcol,'Fontsize',12);
      for elem = 1:nel
          Ycoord = xyY(evtY(elem,:),:);
          plot([Ycoord(2,1) Ycoord(1,1)], [Ycoord(2,2) Ycoord(1,2)],'-.','Color',redcol);    % Line Y_2 - Y_1
          plot([Ycoord(2,1) Ycoord(3,1)], [Ycoord(2,2) Ycoord(3,2)],'-.','Color',redcol);    % Line Y_2 - Y_3    
          plot([Ycoord(3,1) Ycoord(1,1)], [Ycoord(3,2) Ycoord(1,2)],'-.','Color',redcol);    % Line Y_3 - Y_1
      end
      
  elseif reftype == 2
      %
      % Bisec3 subdivision
      %
      % Write vertices' numbers
      text(xyY(:,1),xyY(:,2),vtxnum,'Color',darkGreen,'Fontsize',12);
      for elem = 1:nel
          Xcoord = xy(evt(elem,:),:);
          Ycoord = xyY(evtY(elem,:),:);
          plot([Ycoord(2,1) Ycoord(1,1)], [Ycoord(2,2) Ycoord(1,2)],'-.','Color',darkGreen); % Line Y_2 - Y_1
          plot([Ycoord(2,1) Xcoord(2,1)], [Ycoord(2,2) Xcoord(2,2)],'-.','Color',darkGreen); % Line Y_2 - X_2
          plot([Ycoord(2,1) Ycoord(3,1)], [Ycoord(2,2) Ycoord(3,2)],'-.','Color',darkGreen); % Line Y_2 - Y_3    
      end
  else 
      error('Third input has to be either 1 or 2!');
  end

  axis square;
  set(gca,'XTick',[],'YTick',[],    'Fontsize',32); 
  set(gcf,'units','normalized','Position',[0.25 0.05 0.5 0.7]);
  axis off;
  hold off;

end  % end function