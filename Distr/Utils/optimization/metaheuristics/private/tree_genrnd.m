function tree = tree_genrnd(maxtreedepth,symbols);
        %Generates a random tree
        % tree = tree_genrnd(maxtreedepth,symbols);
        %   tree <- output tree
        %   maxtreedepth -> maximum tree depth
        %   symbols -> cell arrays of operator and terminator node strings
        %

        % (c) Janos Madar, University of Veszprem, 2005

        %Random filling
        nn = [length(symbols{1}), length(symbols{2})];
        n = 2^floor(maxtreedepth)-1;
        vt = zeros(n,1);
        vn = zeros(n,1);

        % Add the PLUS function
        i=1;
        vt(i)=1;
        vn(i)=1;

        % Add a BaseVector from symbols{4}(4) to symbols{4}(5)
        i=2;
        vt(i)=2;
        vn(i)=symbols{4}(4) + randi(symbols{4}(5) - symbols{4}(4) - 1);   %randi(18-9) + 8;

        %for i=1:1,
    %        [vt(i) vn(i)] = [1, 9]; %%tree_genrndsymb(0,nn);
        %end

		not_done = true;
        while (not_done || tree_size(tree) < 5)
			not_done = false;
            for i=3:(n-1)/2,
              [vt(i) vn(i)] = tree_genrndsymb(1/2,nn);
            end
            for i=(n+1)/2:n,
              [vt(i) vn(i)] = tree_genrndsymb(1,nn);
            end
            %Result
            tree.maxsize = n;
            tree.nodetyp = vt;
            tree.node = vn;
            tree.param = zeros(floor((tree.maxsize+1)/2),1);
            tree.paramn = 0;
        end;
end

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
end

