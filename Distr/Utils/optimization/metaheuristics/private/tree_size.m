function [n,vix] = tree_size(tree);
%Gets the size of tree
%  [n,vix] = tree_size(tree)
%    n <- size of tree
%    vix <- index vector of nodes
%    tree -> tree
%

% (c) Janos Madar, University of Veszprem, 2005

%Walking
nt = tree.nodetyp;
v = [1];
i = 1;
n = 1;
while i<=n,
 if nt(v(i))==1,
  v = [v, v(i)*2, v(i)*2+1];
 end
 i = i+1;
 n = length(v);
end
%Result
n = length(v);
vix = v;
