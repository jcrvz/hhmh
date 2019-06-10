%______________________________________________________________________________________________
%  DRONE SQUADRON OPTIMIZATION Algorithm (DSO) toolbox                                                            
%  Source codes demo version 1.0
%                                                                                                     
%  Developed in Octave 3.8. Compatible with MATLAB 2011
%                                                                                                     
%  Author and programmer: Vinicius Veloso de Melo
%                                                                                                     
%         e-Mail: dr.vmelo@gmail.com                                                             
%                 vinicius.melo@unifesp.br                                               
%                                                                                                     
%       Homepage: http://www.sjc.unifesp.br/docente/vmelo/                                                         
%                                                                                                     
%  Main paper:                                                                                        
%  de Melo, V.V. & Banzhaf, W. Neural Comput & Applic (2017). doi:10.1007/s00521-017-2881-3
%  link.springer.com/article/10.1007/s00521-017-2881-3
%_______________________________________________________________________________________________


%%%%%  GPOLS CONFIG

function [config]=InitGPConfig (popsize, maxtreedepth, use_known_perturbation)

  %GP equation symbols

  % Functions with arity 1 (take a single parameter)
  singleParFunctions = { ...
    'psqrt'; ...
    'square'; ...
    %'abs'; ...
    %'neg'; ...
  };

  % Functions with arity 2 (take two parameters)
  doubleParFunctions = {
   'plus'; ...
   'times'; ...
   'minus'; ...
   %'pdiv'; ...
   %'avg'; ...
  };


  constants = {
          'C1'; ...
          'C2'; ...
          'C3'; ...
          'matInterval'; ...      %% range of the problem's bounds, for instance [-10, 10] = 20
  };

  randomNumbers = {
        'U0_1(N)'; ...
        'U05_1(N)'; ...
        'N0_1(N)'; ...
        '(N05_01(N))'; ...
        '(N0_01(N))'; ...
  };

  calculatedWeights = {
        'Sigma(LocalBestPosition, LocalBestOFV, N, matInterval)'; ...
        %'Sigma(LocalBestPosition, LocalBestOFV, N, matInterval)'; ...
        'myeye(Nvars, randi(Nvars, 1, N))'; ...
        'repmat(std(LocalBestPosition), N, 1)'; ...
        'repmat(std(LocalBestPosition(posPBest, :)), N, 1)'; ...
  };

  % Set the possible Departure coordinates
  departureCoordinates = {
          'GlobalBestPosition'; ...
          'CMA (LocalBestPosition, LocalBestOFV, N)'; ...
          'LocalBestPosition'; ...
          'LocalBestPosition (posPBest, : )'; ...
          'LocalBestPosition (a1, : )'; ...     %% doesn't allow repetition
  };

  % Set the other possible coordinates to be used in Offset calculation
  differentialVectors = {
          'CMA (LocalBestPosition, LocalBestOFV, N)'; ...
          'CurrentPosition'; ...
          'Shift'; ...
          'CurrentPosition (a2, : )'; ...    %% allows repetition
          'CurrentPosition (a3, : )'; ...       %% doesn't allow repetition         
          'LocalBestPosition'; ...
          'LocalBestPosition (posPBest, : )'; ...  %% only the pBest positions
          'LocalBestPosition (a2, : )'; ...     %% doesn't allow repetition         
          'LocalBestPosition (a3, : )'; ...     %% doesn't allow repetition
          'GlobalBestPosition'; ...
          'Opposition (LocalBestPosition(a2, :), sumIntervals)'; ...
          'Opposition (CurrentPosition, sumIntervals)'; ...          
  };

  % Symbols can have arity 2 or 1
  symbols{1} = [doubleParFunctions; singleParFunctions ];
  symbols{3} = num2cell([2*ones(1, length(doubleParFunctions)), ones(1, length(singleParFunctions))]);

  symbols{2} = [constants; randomNumbers; calculatedWeights; departureCoordinates; differentialVectors];
  symbols{4} = [1, 1+length(constants), 1+length(constants)+length(randomNumbers), 1+length(constants)+length(randomNumbers)+length(calculatedWeights), ...
      1+length(constants)+length(randomNumbers)+length(calculatedWeights)+length(departureCoordinates)];

      
  config.symbols = symbols;
  config.popsize = popsize;
  config.maxtreedepth = maxtreedepth;
  config.popu = gpols_init(popsize,maxtreedepth,symbols);

  %if (loadGraph)
  %    for (i=1:popsize)
  %        config.popu.chrom{i}.tree = tree_mutate_graph(config.popu.chrom{i}.tree, symbols, graph, max_subtree_depth);
  %    end
  %end

  if (use_known_perturbation)

          %% rand/1
          config.popu.chrom{1}.tree.nodetyp(:) = 2;
          config.popu.chrom{1}.tree.nodetyp(1:3) = [1;2;1];
          config.popu.chrom{1}.tree.nodetyp(6) = 2;
          config.popu.chrom{1}.tree.nodetyp(7) = 1;
          config.popu.chrom{1}.tree.nodetyp(14:15) = 2;

          config.popu.chrom{1}.tree.node(:) = 2;
          config.popu.chrom{1}.tree.node(1:3) = [getSymbolPos(symbols, 'plus', 1); getSymbolPos(symbols, 'LocalBestPosition (a1, : )', 2); getSymbolPos(symbols, 'times', 1)];
          config.popu.chrom{1}.tree.node(6) = getSymbolPos(symbols, 'C1', 2);
          config.popu.chrom{1}.tree.node(7) = getSymbolPos(symbols, 'minus', 1);
          config.popu.chrom{1}.tree.node(14) = getSymbolPos(symbols, 'LocalBestPosition (a2, : )', 2);
          config.popu.chrom{1}.tree.node(15) = getSymbolPos(symbols, 'LocalBestPosition (a3, : )', 2);

      if (popsize > 100)
          %% CMA
          config.popu.chrom{2}.tree.nodetyp(:) = 2;
          config.popu.chrom{2}.tree.nodetyp(1:3) = [1;2;2];
          config.popu.chrom{2}.tree.nodetyp(6) = 2;
          config.popu.chrom{2}.tree.nodetyp(7) = 1;
          config.popu.chrom{2}.tree.nodetyp(14:15) = 2;

          config.popu.chrom{2}.tree.node(:) = 2;
          config.popu.chrom{2}.tree.node(1:3) = [getSymbolPos(symbols, 'plus', 1); getSymbolPos(symbols, 'CMA (LocalBestPosition, LocalBestOFV, N)', 2); getSymbolPos(symbols, 'Sigma(LocalBestPosition, LocalBestOFV, N, matInterval)', 2)];
      end

          % 2X rand/1
          %config.popu.chrom{2} = config.popu.chrom{1};

          %% best/1
          % config.popu.chrom{3}.tree.nodetyp(:) = 2;
          % config.popu.chrom{3}.tree.nodetyp(1:3) = [1;2;1];
          % config.popu.chrom{3}.tree.nodetyp(6) = 2;
          % config.popu.chrom{3}.tree.nodetyp(7) = 1;
          % config.popu.chrom{3}.tree.nodetyp(14:15) = 2;
          % 
          % config.popu.chrom{3}.tree.node(:) = 2;
          % config.popu.chrom{3}.tree.node(1:3) = [getSymbolPos(symbols, 'plus', 1); getSymbolPos(symbols, 'GlobalBestPosition', 2); getSymbolPos(symbols, 'times', 1)];
          % config.popu.chrom{3}.tree.node(6) = getSymbolPos(symbols, 'C1', 2);
          % config.popu.chrom{3}.tree.node(7) = getSymbolPos(symbols, 'minus', 1);
          % config.popu.chrom{3}.tree.node(14:15) = getSymbolPos(symbols, 'LocalBestPosition (randperm(N), : )', 2);
          % 
  end


  config.opt = [0.8 0.0 1.0 2 1 0.2 30 0.05 0 0];

  config.maxgen = 1;
