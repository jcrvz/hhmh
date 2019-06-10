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
% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% ObjectiveFunction = @YourCostFunction
% Nvars = number of your variables
% LB=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% UB=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% setup = DSO configuration
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run DSO: [X,Fx, CurveOFV, NFE, StopMessage, HistoryRanks, Firmwares, LocalBestPosition] = DSO(ObjectiveFunction, Nvars,LB,UB, setup)
%______________________________________________________________________________________________


function [X,Fx, StatsCurvesOFV, NFE, StopMessage, HistoryRanks, Firmwares, LocalBestPosition] = DSO(ObjectiveFunction, Nvars, LB, UB, setup)

	stderr = 1;

        more off

        global Teams_OFV

        addpath ('../GPOLS')

        %% Check Input                   
        if nargin < 4
            fprintf(stderr, '\n!!! Call dso(Nvars,LB,UB) or dso(Nvars,LB,UB, setup)\n');
            error('DSO : Not Enough Input Arguments!')
        end

        if ~isequal(Nvars,size(LB,2),size(UB,2))
            error('DSO : Invalid Arguments : isequal(Nvars,size(LB,2),size(UB,2)) should be true ');
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Default Parameters     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        C1 = 0.8;
        C2 = 0.012;
        C3 = 0.72;
        p=0.25;
        Pacc = 0.1;
        ConvThres = 1e-2;
        MaxIterations = 10;
        maxtreedepth = 10;
        Command_Center_Iter = 100;
        use_known_perturbation = 1;
        N = 30;
        N_Teams = 4;
        ReportLag = 100;
        NFE = 0;
        VTR = -Inf;
        MaxStagnation = 10;
        Toolbox = false;
        StopMessage = 'Stopping DSO : Reached Maximum Number of Iterations';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Update setup
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(setup, 'C1'), C1 = setup.C1 ; end
        if isfield(setup, 'C2'), C2 = setup.C2 ; end
        if isfield(setup, 'C3'), C3 = setup.C3 ; end            
        if isfield(setup, 'N'), N = setup.N ; end
        if isfield(setup, 'N_Teams'), N_Teams = setup.N_Teams ; end
        if isfield(setup, 'MaxIterations'), MaxIterations = setup.MaxIterations ; end
        if isfield(setup, 'VTR'), VTR = setup.VTR ; end
        if isfield(setup, 'ReportLag'), ReportLag = setup.ReportLag ; end
        if isfield(setup, 'MaxStagnation'), MaxStagnation = setup.MaxStagnation ; end          
        if isfield(setup, 'Command_Center_Iter'), Command_Center_Iter = setup.Command_Center_Iter ; end  
        if isfield(setup, 'Pacc'), Pacc = setup.Pacc ; end  
        if isfield(setup, 'ConvThres'), ConvThres = setup.ConvThres ; end          
        if isfield(setup, 'Toolbox'), Toolbox = setup.Toolbox; end          

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Stagnation=0;

        disp ('');
        fprintf('==================================\n', ...
        'Using the following setup:\n ',...
         'C1 = %f\n ',...
         'C2 = %f\n ',...
         'C3 = %f\n ',...       
         'Maximum Iterations = %d\n ',...
         'Squadron size = %d\n ',...
         'Number of teams = %d\n ',...
         'Report Lag = %d\n ',...
        '==================================\n\n', C1, C2, C3, MaxIterations, N, N_Teams, ReportLag);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% AUXILIARY FUNCTIONS %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        U0_1 = @(N) repmat(randraw('unif', [0, 1], N, 1), 1, Nvars);
        U05_1 = @(N) repmat(randraw('unif', [0.5, 1], N, 1), 1, Nvars);
        N0_1 = @(N) repmat(randraw('normal', [0, 1], N, 1), 1, Nvars);
        N05_01 = @(N) repmat(randraw('normal', [0.5, 0.1], N, 1), 1, Nvars);
        N0_01 = @(N) repmat(randraw('normal', [0, 0.1], N, 1), 1, Nvars);

        Opposition = @(LocalBestPosition, sumIntervals) sumIntervals - LocalBestPosition;

        avg = @(a, b) (a+b)/2;        
        %pdiv = @(a,b) (abs(b) > 1e-20) .* a./b;
        plog = @(a) (a>0) .* log(a);
        square = @(a) a .^ 2;
        psqrt = @(a) sign(a) .* sqrt(abs(a));
        neg = @(a) (-1) * a;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        %%%%% GENERATE INITIAL LARGE (N * N_Teams) SAMPLE FROM THE LANDSCAPE
        %% Initial Position
        CurrentPosition = zeros(N*N_Teams,Nvars); % Initial Position
        for i = 1:Nvars
            CurrentPosition(:,i) = randraw('unif', [LB(i),UB(i)],N*N_Teams,1);
        end


        %% Evaluate Initial Position
        CurrentOFV = zeros(N*N_Teams,1); % OFV Value
        for i = 1:N*N_Teams
            CurrentOFV(i) = ObjectiveFunction(CurrentPosition(i,:));
            NFE = NFE + 1;
        end

        %% Select the N best solutions
        [sorted, indices] = sort(CurrentOFV);

        %% Update Local Best
        LocalBestPosition = CurrentPosition( indices(1:N), : ); % Local Best
        LocalBestOFV = CurrentOFV ( indices(1:N) );

        pbest = ceil(N*p);
        posPBest = randi(pbest, 1, N);   %indices(randi(pbest, 1, N));        

        CurrentPosition = LocalBestPosition ;
        CurrentOFV = LocalBestOFV;

        %% Update Global Best
        [GlobalBestOFV, index] = min(LocalBestOFV);
        GlobalBestPosition = repmat(LocalBestPosition(index,:),N,1); % Global Best

        FinalBestOFV = GlobalBestOFV;
        FinalBestSolution = GlobalBestPosition(1,:);
        FinalBestNFE = NFE;


        %% Initial Shift ( Velocity )
        Shift = CurrentPosition .* rand(N,Nvars) ;

        matInterval = repmat(abs(UB-LB),N,1);
        sumIntervals  = repmat(abs(UB+LB),N,1);

        %% Initialize Teams and Firmwares

        Teams_Positions = cell(1, N_Teams); % matrix of solutions
        Teams_Positions_Evaluation = Inf * ones(N, N_Teams);  % matrix of objective function values
        Teams_Ranking = Teams_Positions_Evaluation;  % matrix of objective function values
        Teams_OFV = zeros(1, N_Teams);

        Firmwares = Teams_Positions;

        GPConfig = InitGPConfig(N_Teams, maxtreedepth, use_known_perturbation);

        for i = 1:N_Teams
                Firmwares{i} = tree_stringrc( GPConfig.popu.chrom{i}.tree, 1, GPConfig.symbols );
        end

        fprintf(stderr, 'INITIAL FIRMWARES:');
        disp(Firmwares);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% MAIN LOOP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        StatsCurvesOFV = Inf(MaxIterations, 5);
        
        StatsCurvesOFV(1, :) = [GlobalBestOFV, min(CurrentOFV), mean(CurrentOFV), median(CurrentOFV), max(CurrentOFV)];
        
        %CurveOFV = Inf(1, MaxIterations);
        %CurveOFV(1) = GlobalBestOFV;

        PreviousGlobalBestOFV = Inf;

        rot = (0:1:N-1);               % rotating index array (size NP)
        rotd= (0:1:Nvars-1);                % rotating index array (size D)        

        HistoryRanks = Inf(MaxIterations, N_Teams);


        %% Start search process
        for Iter = 2:MaxIterations

            if abs(FinalBestOFV) <= VTR
                fprintf(stderr, 'VTR FOUND!!! ');
                break;
                continue;
            end

            if NFE > setup.MaxEvaluations
                fprintf(stderr, '!!! NFE > setup.MaxEvaluations !!!\n');
                break;
            end


            ind = randperm(4);              % index pointer array

            a1  = randperm(N);             % shuffle locations of vectors
            rt = rem(rot+ind(1),N);        % rotate indices by ind(1) positions
            a2  = a1(rt+1);                 % rotate vector locations
            rt = rem(rot+ind(2),N);
            a3  = a2(rt+1);                
            rt = rem(rot+ind(3),N);
            a4  = a3(rt+1);               
            rt = rem(rot+ind(4),N);
            a5  = a4(rt+1);              

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% MOVE THE TEAMS USING THE Firmwares
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            violations = zeros(1, N_Teams);
            
            for team = 1:N_Teams
            
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% CALCULATE 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                %% Update Velocity and Position
            		Trial_Position = eval(Firmwares{team});

                mui = rand(N, Nvars) < 0.4 + rand*0.5;          % all random numbers < CR are 1, 0 otherwise

                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% RECOMBINATION
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %% Choose a recombination operator
                opCrossOver = randi(3);

                if (opCrossOver == 1)    % BINOMIAL crossover
                    mpo = mui < 0.5;                % inverse mask to mui
                    jrand = randi(Nvars);
                    New_Team_Position = LocalBestPosition.*mpo + Trial_Position.*mui;
                    
                elseif (opCrossOver == 2)   % EXPONENTIAL CROSSOVER
                    mui=sort(mui');           % transpose, collect 1's in each column
                    for i=1:N
                      n=floor(rand*Nvars);
                      if n > 0
                         rtd = rem(rotd+n,Nvars);
                         mui(:,i) = mui(rtd+1,i); %rotate column i by n
                      end
                    end
                    mui = mui';           % transpose back
                    mpo = mui < 0.5;                % inverse mask to mui                    

                    jrand = randi(Nvars);

                    New_Team_Position = LocalBestPosition.*mpo + Trial_Position.*mui;
                    
                else  % NO CROSSOVER
                    New_Team_Position = Trial_Position;
                end
                   
                   
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% BOUND CORRECTION
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %% Choose a bound correction approach
                opBounds = randi(3);
                
                %% Disable bound correction if it is in the setup
                if (setup.checkBounds == 0)
                    opBounds = -1;
                end

                if (opBounds == 1)
                    for i = 1:Nvars
                        indexes = find(New_Team_Position(:,i) < LB(i));
                        
                        violations(team) = violations(team) + abs(sum(LB(i)-New_Team_Position(indexes,i)));
                        New_Team_Position(indexes,i) = LB(i);
                        indexes = find(New_Team_Position(:,i) > UB(i));
                        
                        violations(team) = violations(team) + abs(sum(New_Team_Position(indexes,i)-UB(i)));
                        New_Team_Position(indexes,i) = UB(i);

                    end
                elseif (opBounds == 2)
                    for i = 1:Nvars
                        indexes = find(New_Team_Position(:,i) < LB(i));
                        violations(team) = violations(team) + abs(sum(LB(i)-New_Team_Position(indexes,i)));
                        New_Team_Position(indexes,i) = LB(i) + (rand(length(indexes),1)*0.1*matInterval(i));
                        indexes = find(New_Team_Position(:,i) > UB(i));
                        violations(team) = violations(team) + abs(sum(New_Team_Position(indexes,i)-UB(i)));
                        New_Team_Position(indexes,i) = UB(i) - (rand(length(indexes),1)*0.1*matInterval(i));
                    end                    
                elseif (opBounds == 3)
                    for i = 1:Nvars
                        indexes = find(New_Team_Position(:,i) < LB(i));
                        violations(team) = violations(team) + abs(sum(LB(i)-New_Team_Position(indexes,i)));
                        New_Team_Position(indexes,i) = LB(i) + abs(rem(New_Team_Position(indexes,i), matInterval(indexes, i)));

                        indexes = find(New_Team_Position(:,i) > UB(i));
                        violations(team) = violations(team) + abs(sum(New_Team_Position(indexes,i)-UB(i)));
                        New_Team_Position(indexes,i) = UB(i) - abs(rem(New_Team_Position(indexes,i), matInterval(indexes, i)));

                    end
                else
                    
                end

                    
                Teams_Positions{team} = New_Team_Position;                

                %% Evaluate New_Team_Position Position
                for idx = 1:N
                    Teams_Positions_Evaluation(idx, team) = ObjectiveFunction(New_Team_Position(idx,:)); % + violations;
                end
                NFE = NFE + N;
                
            end  %  for i = 1:N_Teams

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%  COMMAND CENTER
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Select for next gen the best solutions by row
            for idx = 1:N
                [sorted, indices] = sort(Teams_Positions_Evaluation(idx, :));  % Sort the row position idx

                CurrentOFV (idx) = sorted(1);  % Get the minimum value
                CurrentPosition (idx,:) = Teams_Positions{indices(1)}(idx, :) ;  % Get the solution from Team (indices(1)) that gave the minimum value
                
                Teams_Ranking (idx, :) = indices;
            end
            
            Shift = LocalBestPosition - CurrentPosition;

            %%% Acummulate the historical ranks
            for (team = 1:N_Teams)
              if (std(Teams_Positions_Evaluation(:, team)) == 0)    %%% ignore clones
                Teams_Ranking(:, team) = Inf;
              end;
            end;
 
            Teams_OFV = Teams_OFV + mean(Teams_Ranking, 1) + violations;  % average rank of each team

            HistoryRanks (Iter-1, :) = mean(Teams_Ranking, 1);

            %% Update Local Best
            improved = CurrentOFV < LocalBestOFV;

      	    [sorted2, indices2] = sort(CurrentOFV);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ACCEPT WORSE SOLUTIONS IF SEARCH STAGNATED ?  Pacc% OF CHANCE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (Stagnation > MaxStagnation-1),
              fprintf(stderr, ' Stagnated!!! Selecting non-improved solutions.\n');
              indexes = find( (improved + (rand(N, 1) < Pacc))  > 0) ;
              Stagnation = 0;
            else
              indexes = find( improved ) ;
            end;
            
            LocalBestOFV(indexes) = CurrentOFV(indexes);
            LocalBestPosition(indexes,:) = CurrentPosition(indexes,:);
            
            %% Update Global Best
            [GlobalBestOFVNew,index] = min(LocalBestOFV);
            if GlobalBestOFVNew < GlobalBestOFV
                GlobalBestOFV = GlobalBestOFVNew;
                GlobalBestPosition = repmat(LocalBestPosition(index,:),N,1);
            end    


            meanOFV = mean(Teams_Positions_Evaluation, 1);
            stdOFV = std(Teams_Positions_Evaluation);

            [sortedMeans, indices] = sort(meanOFV);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TIME TO UPDATE THE FIRMWARE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (rem (Iter,  Command_Center_Iter) == 0)

                    %% Update FITNESS in GPLab
                    for i = 1:N_Teams

                      GPConfig.popu.chrom{i}.fitness = Teams_OFV(i);

                      %%% Clear the ranks
                      Teams_OFV(i) = 0;

                    end                  

                    %%% GENERATE A NEW FIRMWARE TO REPLACE THE WORST ONE
                    [GPConfig.popu, dummy] = gpols_mainloop(GPConfig.popu, 0, 0, [], GPConfig.opt);

                    for i = 1:N_Teams
                            Firmwares{i} = tree_stringrc( GPConfig.popu.chrom{i}.tree, 1, GPConfig.symbols );
                    end
                    
        end %% Command_Center_Iters

            %% Print Result
            Avg = mean(CurrentOFV);
            Worst = max(CurrentOFV);
            reference = min(CurrentOFV);            
            
            Restart = 0;
            
            if (Iter==1)||(rem(Iter, ReportLag)==0)
                fprintf(stderr,'Iter: %d | F(x) Best: %.3e | Avg: %.3e | Worst: %.3e | Std-dev: %.3e\n', Iter,GlobalBestOFV, Avg, Worst, std(CurrentOFV));
                
                if (Toolbox)
                    str = [cellstr(get(setup.handles.txt_LOG, 'String')); sprintf('Iter: %d | F(x) Best: %.3e | Avg: %.3e | Worst: %.3e | Std-dev: %.3e', Iter,GlobalBestOFV, Avg, Worst, std(CurrentOFV))];
                    
                    set(setup.handles.txt_LOG, 'String', str);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% PLOT
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    semilogy(StatsCurvesOFV, 'LineWidth', 2)

                    % Turn on the grid
                    grid on
                    
                    legend('Global Best','Current Best','Current Average','Current Median','Current Worst');

                    % Add title and axis labels
                    title(sprintf('Value: %.3e     Evaluations: %d', GlobalBestOFV, NFE))
                    xlabel('Iterations')
                    ylabel('log10 (Objective Function Value)')
                    drawnow                    
                end
            end

            
            if ( std(CurrentOFV)  < ConvThres )   %% Shall we restart due to convergence ?
                Restart = 1;
            elseif ( abs(GlobalBestOFV - PreviousGlobalBestOFV) / abs(PreviousGlobalBestOFV)  < 1e-4 )   %% Needs better improvement
                Stagnation = Stagnation + 1;                
            else
                Stagnation = 0;                
            end

            if (Restart == 1)
                    
                    if Restart == 1
                        fprintf(stderr, 'CONVERGENCE DETECTED AT ITER=%d!!!\n', Iter);
                    else
                        fprintf(stderr, 'STAGNATION DETECTED AT ITER=%d!!!\n', Iter);
                    end
                       
                    Stagnation=0;

                    GPConfig = InitGPConfig(N_Teams, maxtreedepth, use_known_perturbation);

                    for i = 1:N_Teams
                        Firmwares{i} = tree_stringrc( GPConfig.popu.chrom{i}.tree, 1, GPConfig.symbols );
                    end

                   
                    fprintf(stderr, ' ## NEW SAMPLE SIZE=%d\n\n', N);

                    %%%%% GENERATE INITIAL LARGE (N * N_Teams) SAMPLE FROM THE LANDSCAPE
                    %% Initial Position and Velocity
                    len = N*N_Teams;
                    CurrentPosition = zeros(len,Nvars); % Initial Position
                    for i = 1:Nvars
                        CurrentPosition(:,i) = randraw('unif',[LB(i),UB(i)],len,1);
                    end

                    %% Evaluate Initial Position
                    CurrentOFV = zeros(len,1); % OFV Value
                    for i = 1:len
                        CurrentOFV(i) = ObjectiveFunction(CurrentPosition(i,:));
                    end
                    NFE = NFE + len;                                        

                    %% Select the N best solutions
                    [sorted, indices] = sort(CurrentOFV);

                    %% Update Local Best
                    LocalBestPosition = CurrentPosition( indices(1:N), : ); % Local Best
                    LocalBestOFV = CurrentOFV ( indices(1:N) );

                    CurrentPosition = LocalBestPosition ;
                    CurrentOFV = LocalBestOFV;

                    %% Update Global Best
                    [GlobalBestOFVNew, index] = min(LocalBestOFV);
                    GlobalBestOFV = GlobalBestOFVNew;
                    GlobalBestPosition = repmat(LocalBestPosition(index,:),N,1);
                    
                    %%% UPDATE MATINTERVAL
                    matInterval = repmat(abs(UB-LB),N,1);
                    sumIntervals  = repmat(abs(UB+LB),N,1);
                    Shift = CurrentPosition .* rand(N,Nvars) ;
        	    rot = (0:1:N-1);               % rotating index array (size NP)
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %%%% END COMMAND CENTER
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%fflush(stderr) 
            [sorted, indices] = sort(CurrentOFV);
            pbest = ceil(length(indices)*p);
            posPBest = indices(randi(pbest, 1, N));

            if GlobalBestOFV < FinalBestOFV
                FinalBestOFV = GlobalBestOFV;
                FinalBestSolution = GlobalBestPosition(1,:);
                FinalBestNFE = NFE;
            end


            PreviousGlobalBestOFV = GlobalBestOFV;

            % CurveOFV(Iter) = GlobalBestOFV;
            
            StatsCurvesOFV(Iter, :) = [GlobalBestOFV, min(CurrentOFV), mean(CurrentOFV), median(CurrentOFV), max(CurrentOFV)];
            
        end


        fprintf(stderr, 'FINAL FIRMWARES:\n');

        for i = 1:N_Teams
          disp(Firmwares{i});
          fprintf(stderr, 'ts=%d\n', tree_size(GPConfig.popu.chrom{i}.tree));
        end


        if (Toolbox)
            
            str = [cellstr(get(setup.handles.txt_LOG, 'String')); sprintf('Iter: %d | F(x) Best: %.3e | Avg: %.3e | Worst: %.3e | Std-dev: %.3e\n', Iter,GlobalBestOFV, Avg, Worst, std(CurrentOFV))];
            set(setup.handles.txt_LOG, 'String', str);
            
            str = cellstr(get(setup.handles.txt_LOG, 'String'));
            if abs(FinalBestOFV) <= VTR
                str = [str; {'>>> VTR FOUND! <<<'}; {''}];
            end
            str = [str; {'========================'}; {'FINAL FIRMWARES:'}; {'========================'}];

            for i = 1:N_Teams
              str = [str; {Firmwares{i}}];
            end

            str = [str; {' '}; {' '}];

            %% set(setup.handles.txt_LOG, 'String', sprintf('%s\nIteration-%d | Best F(x) = %e | Avg F(x) = %e | Worst F(x) = %e  |  Convergence level=%e\n', get(setup.handles.txt_LOG, 'String'), Iter,GlobalBestOFV, Avg, Worst, abs(GlobalBestOFV - PreviousGlobalBestOFV) / abs(PreviousGlobalBestOFV) ));                    
            set(setup.handles.txt_LOG, 'String', str);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% PLOT
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            semilogy(StatsCurvesOFV, 'LineWidth', 2)

            % Turn on the grid
            grid on
            
            legend('Global Best','Current Best','Current Average','Current Median','Current Worst');

            % Add title and axis labels
            title(sprintf('Value: %.3e     Evaluations: %d', GlobalBestOFV, NFE))
            xlabel('Iterations')
            ylabel('log10 (Objective Function Value)')
            drawnow                    
        end

        %% Output
        X = FinalBestSolution;
        Fx = FinalBestOFV;
        NFE = FinalBestNFE;        

end

