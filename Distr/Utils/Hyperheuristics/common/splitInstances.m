function splitInstances(instanceName, nbFiles, seed)
rng(seed);
fprintf("Splitting instances from path %s ...\n", instanceName);

instanceFiles = dir(instanceName);
numInstances = length(instanceFiles)-2;
wholeData = cell(1,numInstances);
wholeMaxWeight = zeros(1,numInstances);
instanceIDs = randperm(numInstances);

fprintf("\t--- Moving training instances ...\n", instanceName);
for idx = 1 : nbFiles
    fileOriginName = [instanceName instanceFiles( (instanceIDs(idx)) +2).name];
    fileDestinyNameTrain = [instanceName '../SplittedInstances/Train/' instanceFiles( (instanceIDs(idx)) +2).name];
    copyfile(fileOriginName, fileDestinyNameTrain);
end

fprintf("\t--- Moving testing instances ...\n", instanceName);
for idx = 1 + nbFiles : numInstances
    fileOriginName = [instanceName instanceFiles( (instanceIDs(idx)) +2).name];
    fileDestinyNameTest = [instanceName '../SplittedInstances/Test/' instanceFiles( (instanceIDs(idx)) +2).name];
    copyfile(fileOriginName, fileDestinyNameTest);
end

end