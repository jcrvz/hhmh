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

function pos=getSymbolPos(symbols, symbolString, type)
    pos=find(strcmpi(symbols{type}, symbolString));
    if (length(pos) < 1)
        fprintf('Searching for symbol %s.\n', symbolString);
        error('SYMBOL NOT FOUND!!!');
    else
        pos=pos(1);
    end
    
