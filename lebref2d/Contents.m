% LEBMESHREF2D
%
% Files
%   adjustunstructmesh    - renumbering nodes (i.e., edges) for unstructured meshes  
%   bisection             - performs bisections of marked elements/edges
%   crackdomain           - structured crack domain (crack on the left) grid generator
%   detailgrid            - linear detail space Y grid generator
%   edge2nodes            - returns the extremal nodes of the edgenum-th of the mesh 
%   getallmarkelem        - recovers the overall set of marked elements (and edges) (avoiding hanging nodes)
%   lebdemo               - mesh refinement demo
%   lebmeshref            - mesh refinement based on longest edge bisection (LEB) algorithm
%   lshapedomain          - structured L-shaped domain grid generator
%   lshapedomain_unstruct - unstructured L-shaped domain grid generator
%   plotmarkedelem        - plots the marked elements of a given mesh 
%   plotmesh              - plots the mesh given in input
%   plotmeshxandy         - plots the mesh and its uniform (red/bisec3) refinement 
%   squaredomain          - structured square domain grid generator
%   unimeshref            - uniformly refines the triangular mesh using either red or bisec3 refinement
