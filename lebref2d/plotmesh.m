function plotmesh(MESH,titleplot,plotelem,plotvtx)
%PLOTMESH plots the mesh given in input
%  input 
% ------------
%  MESH      : mesh data structure
%  titleplot : (optional) figure's title
%  plotelem  : (optional) write elements's numbers (1/0)
%  plotvtx   : (optional) write vertices's number (1/0)
%
% Default mesh colour: blue
%
% ----------------------------------------------
% EXAMPLE 1. Calling without any optional input:
%    plotmesh(MESH);
%
% EXAMPLE 2. Calling with a title and display vertices
%    plotmesh(MESH,'some-mesh',0,1);
%
% EXAMPLE 3. Calling with all inputs
%    [MESH] = lshapedomain;
%    plotmesh(MESH,'Mesh example',1,1);
% ----------------------------------------------
%
% See also PLOTMESHXANDY, PLOTMARKEDELEM
%
% LEBREF2D function; Copyright (c) L. Rocchi  

  if nargin < 4
      % no vertices' numbers
      plotvtx = 0;         
      if nargin < 3
          % no elements' numbers
          plotelem = 0;     
          if nargin < 2
              % no title
              titleplot = '';
              if nargin < 1
                  error('LEBREF2D: at least one input required!');
              end
          end
      end
  end
  
% Extract elements and nodes
  evt  = MESH.elem;
  xy   = MESH.coord;
  nel  = size(evt,1);   % number of elements
  nvtx = size(xy,1);    % number of nodes  

% -------------------------------------------------------------      
% Adjacency matrix
% -------------------------------------------------------------    
% This is a sparse adjacency matrix for the mesh given by evt matrix
  adj  = sparse(nvtx,nvtx);
  evtt = [evt,evt(:,1)];
  for j = 1:3
      ncol1 = evtt(:,j);
      ncol2 = evtt(:,j+1);
      adj = adj + sparse(ncol1,ncol2,1,nvtx,nvtx);
  end

% -------------------------------------------------------------    
% Plot
% -------------------------------------------------------------    
  figure;
  hold on;
  gplot(adj,xy,'b');
 
% Display elements' numbers  
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
      % Write elements' number
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

  axis square;  
  title(titleplot,'Fontsize',17);
  set(gca,'FontSize',17);
  set(gcf,'units','normalized','Position',[0.25 0.05 0.5 0.7]);
  
% Axis ticks
  minx = min(xy(:,1)); maxx = max(xy(:,1)); midpx = 0.5*(minx+maxx); midplx = 0.5*(minx+midpx); midprx = 0.5*(midpx+maxx);
  miny = min(xy(:,2)); maxy = max(xy(:,2)); midpy = 0.5*(miny+maxy); midply = 0.5*(miny+midpy); midpry = 0.5*(midpy+maxy);
  xticks = [minx midplx midpx midprx maxx];
  yticks = [miny midply midpy midpry maxy];
 
  set(gca,'XTick',xticks,'YTick',yticks,   'Fontsize',22);
  axis off;
  hold off;

end  % end function