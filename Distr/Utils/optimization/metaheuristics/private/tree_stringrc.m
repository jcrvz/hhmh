function s = tree_stringrc(tree,ix,symbols)
%Decodes the tree to string
% s = tree_stringrc(tree,ix,symbols)
%   s <- the output string
%   tree -> the tree
%   ix -> index of strating point (the root = 1)
%   symbols -> cell arrays of operator and terminator node strings
% 
% Remark: It is a recursive function.
%

% (c) Janos Madar, University of Veszprem, 2005


if tree.nodetyp(ix)==1 && ix*2+1<=tree.maxsize,

  sleft = tree_stringrc(tree,ix*2,symbols);
  sright = tree_stringrc(tree,ix*2+1,symbols);
  	if tree.nodetyp(ix)==1 && symbols{3}{tree.node(ix)} == 1
		s = strcat(symbols{tree.nodetyp(ix)}{tree.node(ix)}, ...
	  		'(',sright,')');  	
	else
		s = strcat('(',symbols{tree.nodetyp(ix)}{tree.node(ix)}, '(', sleft,',  ', ...
			    	sright,'))');
	end
else
  s = symbols{tree.nodetyp(ix)}{tree.node(ix)};
end


