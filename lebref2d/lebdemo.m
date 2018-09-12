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
% LEBREF2D function; 12 September 2018
% Copyright (c) 2018 L. Rocchi

  if nargin < 1
      % Square domain (default)
      domain = 1;
  elseif nargin >= 2
      error('LEBREF2D: too many input values!');
  end
  
  close all;

  fprintf('\n++++++++++++++++++++++++++++++++++++++++++++');
  fprintf('\n++           LEBMESHREF2D DEMO            ++');
  fprintf('\n++++++++++++++++++++++++++++++++++++++++++++');
  
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
  nel  = size(MESHX.elem,1);
  nvtx = size(MESHX.coord,1);
  
% Randomly marking elements
  markpar = 0.5;
  fprintf('\nMarking parameter: %.2f',markpar);

% Plot initial mesh 
  titlemesh = ['#elements = ',num2str(nel),', #nodes = ',num2str(nvtx)];
  plotmesh(MESHX,titlemesh,1,1);
  fprintf('\nInitial mesh data: <strong>%d</strong> elements, <strong>%d</strong> nodes\n',nel,nvtx);
  pause(2.5);
  
% Number of total refinements
  totref = 10;

  for nref = 1:totref
      %
      % Set of marked elements
      Mele = randperm(nel,ceil(markpar*nel))';
      %
      % New mesh refinement
      fprintf(num2str(repmat('-',1,44)));
      fprintf('\nMesh refinement %d',nref);
      meshreftime = tic;
      %
      % Mesh refinement
      [MESHX,MMele] = lebmeshref(MESHX,Mele); 
      %
      % Print data
      fprintf(' (%.6f sec)',toc(meshreftime));
      nel  = size(MESHX.elem,1);
      nvtx = size(MESHX.coord,1);
      fprintf('\n-> elements marked:  %.d',length(Mele));
      fprintf('\n-> elements refined: %.d',length(MMele));
      fprintf('\nNew mesh data: <strong>%d</strong> elements, <strong>%d</strong> nodes\n',nel,nvtx);
      %
      % Plot refined mesh
      titlemesh = ['#elements = ',num2str(nel),', #nodes = ',num2str(nvtx)];
      plotmesh(MESHX,titlemesh);
      pause(1.5);
  end
  fprintf('++++++++++++++++++++++++++++++++++++++++++++\n');
  fprintf('++              END OF DEMO               ++\n');
  fprintf('++++++++++++++++++++++++++++++++++++++++++++\n\n');
 
end % end function