function popu = gpols_init(popusize,maxtreedepth,symbols);
%Initializes population variable
% popu = gpols_init(popusize,maxtreedepth,symbols)
%   popu <- generated individuals (population variable)
%   popusize -> number of individulas (size of population)
%   maxtreedepth -> maximum tree depth
%   symbols -> cell arrays of operator and terminator node strings
%
% E.g.
%   symbols{1} = {'+','*');
%   symbols{2} = {'x1','x2','x3'};
%   popu = gpols_init(20,5,symbols);
%
% Remark: This function does not evaluate the individuals
%

% (c) Janos Madar, University of Veszprem, 2005

popu.generation = 1;
popu.symbols = symbols;
popu.size = popusize;
for i = 1:popusize,
  popu.chrom{i}.fitness = 0;
  popu.chrom{i}.mse = 0;
  popu.chrom{i}.tree = tree_genrnd(maxtreedepth,symbols);
end
