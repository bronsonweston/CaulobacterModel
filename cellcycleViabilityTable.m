function [Viable,G1Arrest,G2Arrest,UnknownArrest] = cellcycleViabilityTable(celltype,CtrAExpProf,k,paramChangesStand,paramChangesMult,Name)
%Returns the fraction of cells that are viable or are arrested in G1, G2 or
%unknown. G1 arrest=arrested with 1 chromosome in cell. G2 arrest= arrested
%with 2 chromosomes in cell. Unknown arrest= arrested with any other number
%of chromosomes.
allMutants = {'WT','CtrAdelta3omega','CtrAD51E','CtrAD51Edelta3omega','SM921','deltaGcrA','deltaccrM', ...
    'LS1','PpleC::Tn','divLts','divL(A601L)','divL(Y550F)', 'PpleC::Tn & divL(Y550F)', ...
    'deltadivJ','PpleC::Tn & deltadivJ','deltaSpmX','PpleC::Tn & deltaSpmX', 'PdivK::Tn', ...
    'divKcs','divKxyl','cdG0','PdivK::Tn & cdG0','deltapopA','deltaPleD','deltapopA & PdivK::Tn' ...
    'deltapopA & deltaPleD','deltadgcB','deltaPdeA','deltapopA & deltaPleD & PdivK::Tn', ...
    'cckA(Y514D)','cckA(Y514D) & PdivK::Tn','deltaPodJ','deltahdaA','deltaPleD & PdivK::Tn', 'deltaRcdA'};
myFolder=strcat('ParamCatalog/',CtrAExpProf); 
filePattern = fullfile(myFolder, '*.mat');
matFiles = dir(filePattern);
CountStorage=zeros(length(allMutants),4);
Samples=randsample(length(matFiles),k);
parfor m=1:length(allMutants)
    strain=allMutants{m};
    G1Arrest=0;
    G2Arrest=0;
    UnkownArrest=0;
    Viable=0;
    for i=length(Samples):-1:1
%         global checkmeouty0
        tic
        file=Samples(i);
        baseFileName = matFiles(file).name;
        fullFileName = fullfile(myFolder, baseFileName);
        fprintf(1, 'Now reading %s\n', fullFileName);
        matData = load(fullFileName);
        y0=matData.y0;
        try
            Pset=matData.Pbest;
        catch
            Pset=matData.Pset;
        end
%         [tout, teout, ieout, yout, ~] = CauloSim(celltype,strain, 1600,Pset,y0.');
%         y0=checkmeouty0;
        try
            Pset(paramChangesStand(:,1))=paramChangesStand(:,2);
        catch
        end
        try
            Pset(paramChangesMult(:,1))=paramChangesMult(:,2).'.*Pset(paramChangesMult(:,1));
        catch
        end
        try
            [tout, teout, ieout, yout, ~] = CauloSim(celltype,strain, 1600,Pset,y0.');
            [~,y,cellarrest] = getLastCycle(tout,yout,teout,ieout);
        catch
            UnkownArrest=UnkownArrest+1;
            continue
        end
        if cellarrest == false
            Viable=Viable+1;
        elseif y(end,49)==1
            G1Arrest=G1Arrest+1;
        elseif y(end,49)==2
            G2Arrest=G2Arrest+1;
        else
            UnkownArrest=UnkownArrest+1;
        end
    end
    CountStorage(m,:)=[Viable,G1Arrest,G2Arrest,UnkownArrest];
end

%         nameofFile=strcat(Name,'_',CtrAExpProf,'_',celltype);
%         nameofFile=strcat(nameofFile,'.mat');
%     TotalEvals=k-counter;
%     ArrestFraction=counter/(length(tdivDifList)+counter);
%     WTarrestFraction=(k-(length(tdivDifList)+counter))/k;
RowNames=allMutants;
VariableNames={'Strain','Viable','G1Arrest','G2Arrest','UnknownArrest'};
Viable=CountStorage(:,1)/k; G1Arrest=CountStorage(:,2)/k;G2Arrest= CountStorage(:,3)/k;UnknownArrest= CountStorage(:,4)/k;
T = table(RowNames.',Viable,G1Arrest,G2Arrest,UnknownArrest,...
    'VariableNames',VariableNames)
fullFileName = ['Figures/', Name,'.xlsx'];
writetable(T,fullFileName);
end



