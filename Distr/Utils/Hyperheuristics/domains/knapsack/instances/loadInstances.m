function [wholeData, wholeMaxWeight] = loadInstances(instanceName)
% Input: instanceName
% If instanceName is a folder, Then the script loads all files in folder
% Else the script loads the file
%
% File has the extension .kp and is CSV formatted, where the first line
% contains the wholeMaxData. First and second columns are the weight and
% profit, respectively.
%
% Output: wholeData <- cell ( 1 x numOfInstances )
% Each element corresponds to each instance dataset
%
% Output: wholeMaxWeight <- array ( 1 x numOfInstances )
% Each element corresponds to the first element of the second column  of
% each instance dataset: the max weight

dirMode = false;

fprintf("Loading data with instanceMode = %s ...\n", instanceName);

if instanceName(end) == '/', dirMode = true; end

if dirMode    
    instanceFiles = dir(instanceName);
    numInstances = length(instanceFiles)-2;
    wholeData = cell(1,numInstances);
    wholeMaxWeight = zeros(1,numInstances);
    for idx = 1 : numInstances
        data = csvread([instanceName instanceFiles(idx+2).name]);
        items               = data(2:end, :); % (profit, weight)
        wholeData{idx}      = items;
        wholeMaxWeight(idx) = data(1,2);
    end
else
    data = csvread(instanceName);
    items           = data(2:end, :); % (profit, weight)
    wholeData       = items;
    wholeMaxWeight  = data(1,2);
end
end