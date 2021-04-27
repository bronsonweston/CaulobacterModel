function [newMap] = establishNewMap(MapName,TargetFolder,ParamSetChange,k)
MapName=strrep(MapName,' & ','_');
MapName=strrep(MapName,'&','_');
MapName=strrep(MapName,' + ','_');
MapName=strrep(MapName,'+','_');
filename=['ParamCatalog','/','Maps','.mat'];
Maps=load(filename);

FolderName=strcat('ParamCatalog/',TargetFolder);
filePattern = fullfile(FolderName, '*.mat');
matFiles = dir(filePattern);
% Samples=randsample(length(matFiles),k);
count=0;

MapNames=fieldnames(Maps); map=Maps.(MapNames{1});
newMap = containers.Map('KeyType','char','ValueType','any');
for key=keys(map)
    if strcmp(char(key),'ParamChange')
        newMap(char(key))=ParamSetChange;
    else
        newMap(char(key))=0;
    end
end

Samples=randsample(length(matFiles),k);
for i=1:length(Samples)
    tic
    baseFileName = matFiles(Samples(i)).name;
    fullFileName = fullfile(FolderName, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    matData = load(fullFileName);
    try
        Pset=matData.Pbest;
    catch
        Pset=matData.Pset;
    end
    y0=matData.y0;
    Sbest=matData.Sbest;
    
    try
    Pset(ParamSetChange(:,1))=ParamSetChange(:,2);
    catch
    end
    [Failures] = EvaluateMutants(Pset,y0,1);
    disp( keys(Failures))
    for key=keys(Failures)
        newMap(char(key))=newMap(char(key))+1;
    end
    newMap('ParameterCount')= newMap('ParameterCount')+1;
%     disp( newMap('ParameterCount'))
    toc
    i
end
filename=['ParamCatalog','/','Maps','.mat'];
Maps=load(filename);
Maps.(MapName)=newMap;
save('ParamCatalog/Maps.mat','-struct','Maps')


end

