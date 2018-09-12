function [MMele,MMedge] = getallmarkelem(Mele,evtY,edgelep)
%GETALLMARKELEM recovers the overall set of marked elements (and edges) (avoiding hanging nodes)
%  output 
% ----------
%  MMele   : overall set of marked elements that have to be refined
%  MMedge  : overall set of associated marked edges 
%
%  input 
% ----------
%  Mele    : set of marked elements
%  evtY    : MESHY.elem
%  edgelep : midpoint element-position matrix
%
% Given in input the set of marked elements, the function returns the *overall* 
% set of marked elements/edges so that the their bisection does not produce 
% any hanging node during the mesh refinement 
%
% NOTE that *no mesh* refinement is performed here!
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
% LEBREF2D function; 12 September 2018
% Copyright (c) 2018 L. Rocchi

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
  MMedge = unique( evtY(MMele,2),'stable' );
  
end % end function