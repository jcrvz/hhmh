function [popu,evnum] = gpols_mainloop(popuin,X,Y,Q,opt);
%Run one evolutionary loop, makes the next generation
% [popu,evnum] = gpols_mainloop(popuin,X,Y,Q,opt)
%   popu <- next generation of the population
%   evnum <- number of fun.evaulation (usually number of new individuals)
%   popuin -> the population
%   opt -> options vector, GPOLS-parameters
%   X,Y,Q -> input, output and weighting matrices (see gpols_evaluate)
%
% Remark:
%     opt(1): ggap, generation gap (0-1)
%     opt(2): pc, probability of crossover (0-1)
%     opt(3): pm, probability of mutation (0-1)
%     opt(4): selection type (integer, see gpols_selection)
%     opt(5): rmode, mode of tree-recombination (1 or 2)
%     opt(6): a1, first penalty parameter
%     opt(7): a2, second penalty parameter (0 if there is not penalty)
%     opt(8): OLS treshhold real 0-1 or integer >= 2
%     opt(9): if == 1 -> polynomial evaluation
%     opt(10): if == 1 -> evaluate all indv.s not only new offsprings
%

% (c) Janos Madar, University of Veszprem, 2005

popun = popuin.size;
ggap = opt(1);
pc = opt(2);
pm = opt(3);
tsels = opt(4);
rmode = opt(5);


popu = popuin;

%% Selection
popu.chrom{1}.updated = false;

best = popu.chrom{1}.fitness;
worst = popu.chrom{2}.fitness;  % at least the first team remains unchanged
bestix = 1;
worstix = 2;
for i = 2:popu.size,
  popu.chrom{i}.updated = false;
  if popu.chrom{i}.fitness < best,
    best = popu.chrom{i}.fitness;
    bestix = i;
  end
  if popu.chrom{i}.fitness > worst,
    worst = popu.chrom{i}.fitness;
    worstix = i;
  end  
end


%% New generation
%% Replace only the worst solution with a variant of the best solution
popu.chrom{worstix} = popuin.chrom{bestix};
popu.chrom{worstix}.fitness = Inf;

%mutate tree
tree1 = popu.chrom{bestix}.tree;
tree1 = tree_mutate(tree1,popu.symbols);
popu.chrom{worstix}.tree = tree1;

popu.chrom{worstix}.updated = true;

newix=worstix;

popu.generation = popu.generation+1;
evnum = 1;

return


