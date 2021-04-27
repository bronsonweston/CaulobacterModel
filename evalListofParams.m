function [] = evalListofParams(k)
global checkmeouty0;
%
mutantList = {'WT','CtrAdelta3omega','CtrAD51E','CtrAD51Edelta3omega','SM921','deltaGcrA','deltaccrM', ...
    'LS1','PpleC::Tn','divLts','divL(A601L)','divL(Y550F)', 'PpleC::Tn & divL(Y550F)', ...
    'deltadivJ','PpleC::Tn & deltadivJ','deltaSpmX','PpleC::Tn & deltaSpmX', 'PdivK::Tn', ...
    'divKcs','divKxyl','cdG0','PdivK::Tn & cdG0','deltapopA','deltaPleD','deltapopA & PdivK::Tn' ...
    'deltapopA & deltaPleD','deltadgcB','deltaPdeA','deltapopA & deltaPleD & PdivK::Tn', ...
    'cckA(Y514D)','cckA(Y514D) & PdivK::Tn','deltaPodJ','deltahdaA','deltaPleD & PdivK::Tn', 'deltaRcdA'};
try
    filename=['ParamCatalog','/','Maps','.mat'];
    Maps=load(filename);
catch
    SlowMap = containers.Map('KeyType','char','ValueType','any');
    SlowMap('ParameterCount')=0;
    SlowMap('ParamChange')=[];
%     PartialMap = containers.Map('KeyType','char','ValueType','any');
%     PartialMap('ParameterCount')=0;
%     PartialMap('ParamChange')=[];
    CompleteMap = containers.Map('KeyType','char','ValueType','any');
    CompleteMap('ParameterCount')=0;
    CompleteMap('ParamChange')=[];
    ctrAcoriMap = containers.Map('KeyType','char','ValueType','any');
    ctrAcoriMap('ParameterCount')=0;
    ctrAcoriMap('ParamChange')=[];
    %
    for mut=1:length(mutantList)
        SlowMap(strcat(char(mutantList(mut)),'_SW'))=0; SlowMap(strcat(char(mutantList(mut)),'_ST'))=0;
%         PartialMap(strcat(char(mutantList(mut)),'_SW'))=0; PartialMap(strcat(char(mutantList(mut)),'_ST'))=0;
        CompleteMap(strcat(char(mutantList(mut)),'_SW'))=0; CompleteMap(strcat(char(mutantList(mut)),'_ST'))=0;
        ctrAcoriMap(strcat(char(mutantList(mut)),'_SW'))=0; ctrAcoriMap(strcat(char(mutantList(mut)),'_ST'))=0;
    end
    Maps.('SlowMap')=SlowMap;
%     Maps.('PartialMap')=PartialMap;
    Maps.('CompleteMap')=CompleteMap;
    Maps.('ctrAcoriMap')=ctrAcoriMap;
end

%myFolder='ParameterSetsRaw';
for m=1:3
    if m == 1
        curMap='SlowMap';
        directory= 'ParamCatalog/SlowParams_All';
%     elseif m == 2
%         curMap='PartialMap';
%         directory= 'ParamCatalog/PartialParams_All';
    elseif m == 2
        curMap='CompleteMap';
        directory= 'ParamCatalog/CompleteParams_All';
    elseif m== 3
        curMap='ctrAcoriMap';
        directory= 'ParamCatalog/CtrABindingParams_All';
    end
    filePattern = fullfile(directory, '*.mat');
    matFiles = dir(filePattern);
    Samples=randsample(length(matFiles),k);
    map=Maps.(curMap)
    for i=1:k
        tic
        file=Samples(i);
        baseFileName = matFiles(file).name;
        fullFileName = fullfile(directory, baseFileName);
%         fprintf(1, 'Now reading %s\n', fullFileName);
        matData = load(fullFileName);
        try
            Pset=matData.Pbest;
        catch
            Pset=matData.Pset;
        end
        y0=matData.y0;
        Sbest=matData.Sbest;
%         if Sbest>450
%             delete(fullFileName)
%             continue
%         end
%             paramspec = map('ParamChange');
%             Psetnow=Pset;
%             if ~isempty(paramspec)
%                 Psetnow(paramspec(1))=paramspec(2);
%             end
            %         if strcmp(char(MapNames{m}),'NormalMap')
            %             [Failures,M] = EvaluateMutants(Psetnow,y0,1);
            %         else
            %             [Failures,~] = EvaluateMutants(Psetnow,y0,1);
            %         end
            [Failures] = EvaluateMutants(Pset,y0,1);
            for key=keys(Failures)
                map(char(key))=map(char(key))+1;
            end
            map('ParameterCount')= map('ParameterCount')+1;
            Maps.(curMap)=map;
            
            disp(num2str(((m-1)*k+i)/(3*k)));
            save('ParamCatalog/Maps.mat','-struct','Maps')
            toc
    end
end

end