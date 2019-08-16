function [MMele,MMedge] = getallmarkelem(Mset,evtY,edgelep,imark)
%GETALLMARKELEM recovers the overall set of marked elements (and edges) (avoiding hanging nodes)
%  output 
% ----------
%  MMele   : set of overall marked elements that have to be refined
%  MMedge  : set of overall associated marked edges 
%
%  input 
% ----------
%  Mset    : set of marked elements or edges
%  evtY    : MESHY.elem
%  edgelep : midpoint element-position matrix
%  imark   : 1 or 2 if Mset is a vector of elements or edges, respectively
%
% Given in input the set of marked elements, the function returns the set
% of *overall* marked elements/edges so that the their bisection does not 
% produce any hanging nodes during the mesh refinement.
%
% NOTE that *no* mesh-refinement is performed here!
%
% The function works according to the following procedure:
%
% 1 - Find the set of neighbouring elements sharing the midpoints of the
%     longest edges of the marked elements;
%     These elements are in found in edgelep( evtY(Mele,2) ,[1,2]);
%
% 2 - If this set is exactly Mele, it means that no hanging nodes will be
%     produced and no further elements need to be marked. In particular, 
%     this happens if Mele consists of elements sharing only the longest edge;
%
% 3 - Otherwise, there is some edge's midpoint belonging to an element
%     which has not been marked. Hence, that midpoint will be an hanging node. 
%     In this case, we add the associated element to the vector of marked 
%     elements.
%
% The procedure continues until condition 2 is satisfied.
%
% LEBREF2D function; Copyright (c) L. Rocchi  

  if ~ismember(imark,[1,2])
      error('The fourth argument has to be either equal to 1 or 2!'); 
  end
  
  if imark == 1
      % The input argument is a set of marked elements:
      % - assign Mset to Mele;
      % - set Medge=empty;
      Mele  = Mset;
      Medge = [];      
  else%imark==2
      % The input argument is a set of marked edges:
      % - assign Mset to Medge; 
      % - recover the associated set of "marked elements" which share the marked edges
      Medge = Mset;
      Mele  = edgelep(Medge,[1 2]);
      Mele  = unique( reshape(Mele,2*size(Mele,1),1) );
  end

% Marked elements are saved in a cell: this improve efficiency when new marked 
% elements are added during the loop below since the vector is directly replaced 
% by the new one; in particular, neither preallocation nor 'sparse' is required 
  MMele    = cell(1,1);
  MMele{1} = Mele; 
    
  while true 
      
      % Midpoints of longest edges of the current set of marked elements  
      midptle = evtY(MMele{1},2);  
      
      % Neighbouring elements sharing midptle
      neighElem = edgelep( midptle,[1,2] );
      
      % New marked elements: they are those elements that do not belong to Mele
      % but some of their edges has to be bisected in order to keep conformity
      newmarkelem = unique( neighElem( ~ismember(neighElem , MMele{1}) ) );
      
      if isempty(newmarkelem)
          % The conformity is recovered
          break; 
      end
            
      % Update the current set of marked elements
      MMele{1} = [MMele{1}; newmarkelem];
      
  end
  
% Overall set of marked elements 
  MMele = MMele{1};

% Overall set of marked edges (avoid repetitions)
  MMedge = unique( [Medge; evtY(MMele,2)], 'stable' );
  
end % end function