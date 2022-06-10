function [evt] = bisection(MMele,MMedge,markedge,evt,evtY)
%BISECTION performs bisections of marked elements/edges
%  output 
% -----------
%  evt      : MESHX.elem (after refinement)
%
%  input 
% -----------
%  MMele    : overall set of marked elements that have to be refined
%  MMedge   : overall set of associated marked edges 
%  markedge : marked edge-midpoint-position vector
%  evt      : MESHX.elem (before refinement)
%  evtY     : MESHY.elem
%
% LEBREF2D function; Copyright (c) L. Rocchi  

% Extract marked elements according to the refinement type 1/2/3
  refxelem = sum( ismember( evtY(MMele,:) , MMedge) , 2);
  Mone     = sort( MMele( refxelem == 1 ) );                                     % elements marked for 1 bisection
  Mtwo     = MMele( refxelem == 2 );                                             % elements marked for 2 bisections
  Mthree   = sort( MMele( refxelem == 3 ) );                                     % elements marked for 3 bisections
        
% -----------------------------------------------------------------------------
% STEP 1: refine elements marked for 1 bisection (green refinement)
% -----------------------------------------------------------------------------

  vtx      = evt(Mone,:);                                                         % nodes of elements in Mone
  newnodes = full( markedge( evtY(Mone,2) ) );                                    % midpoints of longest edges of Mone

% Refine elements
  firstChild  = [Mone, vtx(:,2), newnodes, vtx(:,1)];                             % top child elements
  secondChild = [Mone, vtx(:,3), newnodes, vtx(:,2)];                             % bottom child elements

% Save elements
  newrows_one = [firstChild, secondChild];       
  newrows_one = reshape(newrows_one',4,2*length(Mone))';

% Transform the new rows per each child in a cell array: it will be length(Mone) 
% long and each cell is a 2-by-4 matrix (where 2 is due to the 2 children)
  cell_one = mat2cell(newrows_one, 2*ones(1,size(newrows_one,1)/2), 4);
  
  
% -----------------------------------------------------------------------------
% STEP 2: refine elements marked for 2 bisections (blue refinement)
% -----------------------------------------------------------------------------

% Separate the elements marked for blue-left from those marked for blue-right 
% refinement. NOTE that: 
% - blue-left refinement means that (also) the 1st element's edge is bisected
% - blue-right refinement means that (also) the 3rd element's edge is bisected
  checkleft    = repmat( [1 1 0] , [size(Mtwo),1] );
  checkright   = repmat( [0 1 1] , [size(Mtwo),1] );
  markEdgesxel = ismember( evtY(Mtwo,:) , MMedge);
  Mtwoleft     = Mtwo( sum( abs(markEdgesxel - checkleft)  , 2) == 0 );           % marked elements for blue-left refinement
  Mtworight    = Mtwo( sum( abs(markEdgesxel - checkright) , 2) == 0 );           % marked elements for blue-right refinement
  Mtwoleft     = sort( Mtwoleft  );
  Mtworight    = sort( Mtworight );
  
% Allocate an empty vector for all children space
  [firstChildLeft,  secondChildLeft,  thirdChildLeft ,...
   firstChildRight, secondChildRight, thirdChildRight] = deal([]); 
  
  if ~isempty(Mtwoleft) 
      vtxLeft     = evt(Mtwoleft,:);                                              % nodes of elements in Mtwoleft 
      leNodeLeft  = full( markedge( evtY(Mtwoleft,2) ) );                         % midpoints of longest edges of Mtwoleft's elements
      mrkNodeLeft = full( markedge( evtY(Mtwoleft,1) ) );                         % midpoints of 1st edges of Mtwoleft's elements
      % Left-Blue refinement
      firstChildLeft  = [Mtwoleft, vtxLeft(:,2),  leNodeLeft, vtxLeft(:,1)];      % top child elements
      secondChildLeft = [Mtwoleft, vtxLeft(:,2), mrkNodeLeft,   leNodeLeft];      % middle child elements
      thirdChildLeft  = [Mtwoleft,   leNodeLeft, mrkNodeLeft, vtxLeft(:,3)];      % bottom child elements
  end
  
  if ~isempty(Mtworight)
      vtxRight     = evt(Mtworight,:);                                            % nodes of elements in Mtworight
      leNodeRight  = full( markedge( evtY(Mtworight,2) ) );                       % midpoints of longest edges of Mtworigth's elements
      mrkNodeRight = full( markedge( evtY(Mtworight,3) ) );                       % midpoints of 3rd edges of Mtworight's elements
      % Right-Blue refinement
      firstChildRight  = [Mtworight, vtxRight(:,1), mrkNodeRight,   leNodeRight]; % top child elements
      secondChildRight = [Mtworight,   leNodeRight, mrkNodeRight, vtxRight(:,2)]; % middle child elements
      thirdChildRight  = [Mtworight, vtxRight(:,3),  leNodeRight, vtxRight(:,2)]; % bottom child elements
  end   
  
% Save elements  
  newrows_two_left  = [firstChildLeft,  secondChildLeft,  thirdChildLeft];
  newrows_two_left  = reshape(newrows_two_left',4,3*length(Mtwoleft))';
%
  newrows_two_right = [firstChildRight, secondChildRight, thirdChildRight];
  newrows_two_right = reshape(newrows_two_right',4,3*length(Mtworight))';
        
% Transform the new rows per each child in cells array: they will be length(Mtwoleft) 
% and length(Mtworight) long and each cell would be a 3-by-4 matrix 
% (where 3 is due to the 3 children)
  cell_two_left  = mat2cell(newrows_two_left,  3*ones(1,size(newrows_two_left,1)/3),  4);
  cell_two_right = mat2cell(newrows_two_right, 3*ones(1,size(newrows_two_right,1)/3), 4);
  

% -----------------------------------------------------------------------------             
% STEP 3: refine elements marked for 3 bisections (bisec3 refinement)
% -----------------------------------------------------------------------------

  vtx           = evt(Mthree,:);                                                  % nodes of elements in Mthree
  leNode        = full( markedge( evtY(Mthree,2) ) );                             % midpoints of longest edges of Mthree
  mrkNodeLeft   = full( markedge( evtY(Mthree,1) ) );                             % midpoints of 1st edges of Mthree
  mrkNodeRight  = full( markedge( evtY(Mthree,3) ) );                             % midpoints of 3rd edges of Mthree
  
% Refine elements
  firstChild    = [Mthree, vtx(:,1), mrkNodeRight,   leNode];                     % top child elements
  secondChild   = [Mthree,   leNode, mrkNodeRight, vtx(:,2)];                     % first middle child elements
  thirdChild    = [Mthree, vtx(:,2),  mrkNodeLeft,   leNode];                     % second middle child elements
  fourtChild    = [Mthree,   leNode,  mrkNodeLeft, vtx(:,3)];                     % bottom child elements

% Save elements
  newrows_three = [firstChild, secondChild, thirdChild, fourtChild];
  newrows_three = reshape(newrows_three',4,4*length(Mthree))';
  
% Transform the new rows per each child in a cell array: it will be length(Mthree) 
% long and each cell would be a 4-by-4 matrix (where 4 is due to the 4 children)
  cell_three = mat2cell(newrows_three, 4*ones(1,size(newrows_three,1)/4), 4);
  

% -----------------------------------------------------------------------------
% STEP 4: create the new element mapping matrix
% -----------------------------------------------------------------------------    

% New elements are inserted in the current evt in the following way.
% Suppose the 23rd element was marked for 2 refinements. This will produce 3 
% children. Then, the 23rd row in evt should be deleted and 3 new rows 
% corresponding to the 3 children will take its place. The 3 children will 
% stay on the 23rd, 24th, and 25th row. 
% Accordingly, the old 24th row will become the 26th row.
%
% The insertion explained above is done by transforming the evt matrix in a cell array

% Create a cell containing the rows of evt (element)
  evt = [(1:size(evt,1))', evt];
  evtcell = mat2cell(evt, ones(1,size(evt,1)), 4);

% Update elements with the new rows of all children that have to be inserted
  evtcell( Mone      ) = cell_one(:);
  evtcell( Mtwoleft  ) = cell_two_left(:);
  evtcell( Mtworight ) = cell_two_right(:);
  evtcell( Mthree    ) = cell_three(:);
 
% Convert the cell to matrix
  evt = cell2mat( evtcell );  
  evt = evt(:,[2 3 4]);

end % end function