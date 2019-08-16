function lebdemo(domain)
%LEBDEMO mesh refinement demo
%  input
% ---------
%  domain : domain generated
%           1 - square [0,1]^2
%           2 - L-shaped [-1,1]^2 \ (-1,0)^2 (structured)
%           3 - L-shaped [-1,1]^2 \ (-1,0)^2 (unstructured)
%           4 - crack [-1,1]^2 \ (-1,0)x{0} 
%
% Default domain=1 (square)
%
% LEBREF2D function; Copyright (c) L. Rocchi  

  if nargin < 1
      % Square domain (default)
      domain = 1;
  elseif nargin >= 2
      error('LEBREF2D: too many input values!');
  end
  
  close all;

  fprintf('\n++++++++++++++++++++++++++++++++++++++++++++');
  fprintf('\n++            LEBREF2D DEMO               ++');
  fprintf('\n++++++++++++++++++++++++++++++++++++++++++++');
  
% ----------------------------------------------------------  
% Domain generation  
% ----------------------------------------------------------
  if domain == 1
      % Square domain
      fprintf('\n<strong>Square domain</strong>');
      [MESHX] = squaredomain;
  elseif ismember(domain,[2,3])
      % L-shaped domain
      fprintf('\n<strong>L-shaped domain</strong>');
      if domain==2
          [MESHX] = lshapedomain;
      else
          [MESHX] = lshapedomainunstruct;
      end
  elseif domain == 4
      % Crack domain
      fprintf('\n<strong>Crack domain</strong>');
      [MESHX] = crackdomain;
  end
  
% Number of elements  
  nel  = size(MESHX.elem,1);
  
% Number of nodes
  nvtx = size(MESHX.coord,1);
    
% ----------------------------------------------------------  
% Detail grid 
% ----------------------------------------------------------  
  [MESHY,edgelep] = detailgrid(MESHX);
% Number of edges
  nedge = size(MESHY.coord,1);

% ----------------------------------------------------------  
% Marking elements or edges?
% ----------------------------------------------------------  
% Marking elements or edges (1/2): default elements (1)  
  imark = 2;
  if ismember(imark,[1,2])
      if imark == 1
          neledg = nel;
      elseif imark == 2
          neledg = nedge;
      end
  else 
      fprintf('Parameter imark has to be eqither 1 or 2')
  end
    
% Marking parameter for random marking of elements or edges
  markpar = 0.5;
  fprintf('\nMarking parameter: %.2f',markpar);
  
% ----------------------------------------------------------  
% Plot initial mesh 
% ----------------------------------------------------------  
  titlemesh = ['#elements = ',num2str(nel),  ...
               ', #edges = ',num2str(nedge), ...
               ', #nodes = ',num2str(nvtx)];
  plotmesh(MESHX,titlemesh,1,1);
  fprintf('\nInitial mesh data: <strong>%d</strong> elements, <strong>%d</strong> nodes\n',nel,nvtx);
  pause(2.5);
  
% Number of total refinements
  totref = 8;
  
% ----------------------------------------------------------  
% Adaptive mesh-refinements
% ----------------------------------------------------------  
  for nref = 1:totref
      
      fprintf(num2str(repmat('-',1,44)));
      fprintf('\nMesh refinement %d',nref);
      
      % Set of marked elements or edges 
      Mset = randperm( neledg, ceil(markpar*neledg) )';
      
      % Mesh refinement
      meshreftime = tic;
      [MESHX,MESHY,MMele,MMedge,edgelep] = lebmeshref(MESHX,MESHY,edgelep,Mset,imark);       
      fprintf(' (%.5f sec)',toc(meshreftime));
      
      % Print data
      if imark == 1
          fprintf('\n-> marked elements:  %.d' ,length(Mset));
          fprintf('\n-> refined elements: %.d' ,length(MMele));
          fprintf('\n-> (edges bisected:  %.d)',length(MMedge));
          % Update the parameter neledg
          neledg = size(MESHX.elem,1);
      else
          fprintf('\n * marked edges:      %.d' ,length(Mset));
          fprintf('\n * refined edges:     %.d' ,length(MMedge));
          fprintf('\n * (elements refined: %.d)',length(MMele));
          % Update the parameter neledg
          neledg = size(MESHY.coord,1);
      end
      nel   = size(MESHX.elem,1);   % New number of elements
      nvtx  = size(MESHX.coord,1);  % New number of vertices
      nedge = size(MESHY.coord,1);  % New number of edges      
      fprintf('\nNew mesh data:\n'); 
      fprintf(' * <strong>%d</strong> elements\n',nel);
      fprintf(' * <strong>%d</strong> nodes\n',nvtx);
      fprintf(' * <strong>%d</strong> edges\n',nedge);

      % Plot the refined mesh
      titlemesh = ['#elements = ',num2str(nel),  ...
                   ', #edges = ',num2str(nedge), ...
                   ', #nodes = ',num2str(nvtx)];
      plotmesh(MESHX,titlemesh);          
      pause(1.5);
  end
  fprintf('++++++++++++++++++++++++++++++++++++++++++++\n');
  fprintf('++              END OF DEMO               ++\n');
  fprintf('++++++++++++++++++++++++++++++++++++++++++++\n\n');
 
end % end function