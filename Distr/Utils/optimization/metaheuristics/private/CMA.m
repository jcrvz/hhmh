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

% Input: the best solutions
% Output: the multivariate normal distribution matrix following the distribution of the selected best solutions

function mat = CMA (LocalBestPosition, CurrentOFV, N)

    [sorted, indices] = sort(CurrentOFV);

    bestidx = indices (1: ceil(length(indices)*0.5));

    mu = mean(LocalBestPosition (bestidx, : ), 1);

    sigma = cov( LocalBestPosition(bestidx, : ));

    mat = mvnrnd(mu, sigma, N);
end
