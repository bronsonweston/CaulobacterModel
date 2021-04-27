function [pMean,pSTD,pRelSTD,pMin,pMax,T] = getParamStatistics(myFolder,name)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes her

load('paramDictionary.mat','paramDictionary');
filePattern = fullfile(myFolder, '*.mat');
matFiles = dir(filePattern);
ParamMatrix=zeros(length(matFiles),142);
for i=1:length(matFiles)
    baseFileName = matFiles(i).name;
    fullFileName = fullfile(myFolder, baseFileName);
%     fprintf(1, 'Now reading %s\n', fullFileName);
    matData = load(fullFileName);
    try
        Pset=matData.Pbest;
    catch
        Pset=matData.Pset;
    end
    if length(Pset) < 139
        Pset=[Pset,1,1,1,0.98];
    end
    Pset(108)=Pset(108)*Pset(32)/Pset(33);
    Pset(32)=1;
    Pset(33)=1;
    ParamMatrix(i,:)=Pset;
end
pMean = zeros(1,length(keys(paramDictionary)));
pSTD = zeros(1,length(keys(paramDictionary)));
pRelSTD = zeros(1,length(keys(paramDictionary)));
pMin = zeros(1,length(keys(paramDictionary)));
pMax = zeros(1,length(keys(paramDictionary)));

for param=1:length(ParamMatrix(1,:))
    pMean(param)=mean(ParamMatrix(:,param));
    pSTD(param)=std(ParamMatrix(:,param));
    pRelSTD(param)=pSTD(param)/pMean(param);
    pMin(param)=min(ParamMatrix(:,param));
    pMax(param)=max(ParamMatrix(:,param));
end

% TableCells={[] [] [] [] [] []}
TableCells={};
for k=keys(paramDictionary)
    key=char(k);
    p=paramDictionary(key);
%     key
    if isempty(TableCells)
        TableCells={key pMean(p) pSTD(p) pRelSTD(p) pMin(p) pMax(p) p};
    else
        TableCells(length(TableCells(:,1))+1,:)={key pMean(p) pSTD(p) pRelSTD(p) pMin(p) pMax(p) p};
    end
end
clc
T = cell2table(TableCells,...
    'VariableNames',{'Param' 'Avg' 'Std' 'Norm_Std' 'Min' 'Max' 'ParamNum'});
T = sortrows(T,'ParamNum')

if ~isempty(name)
    writetable(T,['Figures/',name, '_Statistics.xlsx']);
end

end