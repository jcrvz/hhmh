% Script for reading 1D BPP instances 
%--- Inputs
%instanceMode: Text providing the location of the file or dir to load. If it is a dir, the final character must be '\'
%--- Outputs
%wholeData: Cell array with the information of each instance
%wholeMaxWeight: Vector with the capacity of bins for each instance
%nbItems: Vector with the number of items in each instance
function [wholeData, wholeMaxWeight, nbItems] = loadInstancesBPP(instanceMode)
 
  fileBase = ['..\domains\binpacking\instances\' instanceMode];
    
  if instanceMode(end) ~= '\'
    isDir = false; 
    nbInstances = 1;
  else
    isDir = true;
    files = dir(fileBase);
    nbInstances = length(files)-2;
  end
  
  allInstances = cell(nbInstances,1);                        
  for idy = 1 :  nbInstances
    if isDir
      allInstances{idy} = [fileBase files(idy+2).name];        
    else
      allInstances{idy} = fileBase;        
    end
  end
    
  wholeData = cell(nbInstances,1);
  wholeMaxWeight = zeros(1, nbInstances);
  nbItems = zeros(1, nbInstances);

  fprintf("Loading data with instanceMode = %s ...\n", instanceMode);

  for idy = 1 : nbInstances    
    data = csvread(allInstances{idy});
    items = data(3:end, :); % (weight)
    wholeData{idy} = items;
    wholeMaxWeight(idy) = data(2,1); % Second row gives the maximum capacity for each bin 
    nbItems(idy) = data(1,1); % First row gives the number of items 
  end
end