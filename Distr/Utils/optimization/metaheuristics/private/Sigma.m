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


% Function modified from the CMA-ES source code

function sigma = Sigma (LocalBestPosition, LocalBestOFV, N, matInterval)
	[sorted, indices] = sort(LocalBestOFV);
    bestidx = indices (1: ceil(length(indices)*.5));
    xmean = mean(LocalBestPosition (bestidx, : ), 1);
    mu = length(bestidx);
    weights = log(mu+0.5)-log(1:mu)'; % muXone array for weighted recombination
	mueff=sum(weights)^2/sum(weights.^2); % variance-effective size of mu    
    sigma = 0.04 * mueff * sqrt(sum(xmean.^2)) / N; % 20D,lam=1000:25e3

    if (sigma > 2) 
    	sigma = 2;
    elseif (sigma < 1e-16)
    	sigma = 1e-16;
    end
    sigma = sigma .* randraw('normal', [0, 1], N, size(LocalBestPosition, 2)) .* matInterval * rand/2;
end
