function tree = tree_mutate(treein,symbols);
%Mutates a tree (mutates one randomly selected node)
% tree = tree_mutate(treein,symbols)
%   tree <- the output tree
%   treein -> the input tree
%   symbols -> cell arrays of operator and terminator node strings
%

% (c) Janos Madar, University of Veszprem, 2005

%Begin

len1 = length(treein.nodetyp);
len2 = length(treein.node);

cont = 0;
not_done=true;
ts = 0; %tree_size(tree);

%% Number of trials to get a different solution
while cont < 3,
	tree = treein;
		nn = [length(symbols{1}), length(symbols{2})];
		%Mutate one node
		[n,v] = tree_size(tree);
		i = v(floor(rand*(length(v))+1));   %% +1 to ignore first PLUS operator
		if i<(tree.maxsize+1)/2 && rand<0.5,
		%if i<(tree.maxsize-1)/2 & rand<0.5,
		  [tree.nodetyp(i) tree.node(i)] = tree_genrndsymb((tree.nodetyp(i)==1),nn);
		else
		  while tree.node(i)==treein.node(i) && tree.nodetyp(i)==treein.nodetyp(i),
		    [tree.nodetyp(i) tree.node(i)] = tree_genrndsymb((tree.nodetyp(i)~=1),nn);
		  end
		end


		%% Repair doubles if needed
		%% for instance: times(C2, C2)  or minus(LocalBestPosition (randperm(N), : ),LocalBestPosition (randperm(N), : ))
	if true,
        	nt = tree.nodetyp;
		v = [1];
		i = 1;
		n = 1;
		while i<=n,
		 if nt(v(i))==1,
		  v = [v, v(i)*2, v(i)*2+1];
		  while (nt(v(i)*2)~=1 && tree.node(v(i)*2) == tree.node(v(i)*2+1))
		  	[dummy tree.node(v(i)*2)] = tree_genrndsymb(1,nn);  %% get a new terminal
		  end
		 end
		 i = i+1;
		 n = length(v);
		end
    	end

		%% Repair Base Vector if needed
		%% 2 = Terminal     symbols{4}(4) = Index of the first Base Vector
		tree.nodetyp(1)=1;
		tree.node(1)=1;

		if tree.nodetyp(2) ~=2 || tree.node(2) < symbols{4}(4) || tree.node(3) < symbols{4}(4) || tree.node(2) >= symbols{4}(5)
		    tree.nodetyp(2)=2;
		    tree.node(2)=symbols{4}(4) + randi(symbols{4}(5) - symbols{4}(4) - 1);
		end

		cont=cont+1;
		ts = tree_size(tree);
	%end


	if ( sum(tree.nodetyp == treein.nodetyp) ~= len1 && sum(tree.node == treein.node) ~= len2 && (ts > 3 && ts < 25))
		break;
	end;

end  %while



%------------------------------------------------------------------
function [nodetyp,node] = tree_genrndsymb(p0,nn)
%Generate a random symbol (tarminate or operate)
%  [nodetyp,node] = tree_genrndsymb(p0,nn)
%   nodetyp,node <- results
%   p0 -> probability of terminate node
%   nn -> vector [number of operators, variables]
%

if rand<p0,
  nodetyp = 2;
else
  nodetyp = 1;
end
node = floor(nn(nodetyp)*rand)+1;
