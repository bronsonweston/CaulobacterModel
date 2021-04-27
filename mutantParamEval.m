function [] = mutantParamEval(celltype,CtrAExpProf,k,paramChangesStand,paramChangesMult,Name)
%   Detailed explanation goes here
%     myFolder='ParameterSetsAnalyzed'; filePattern = fullfile(myFolder, '*.mat');
    myFolder=strcat('ParamCatalog/',CtrAExpProf); filePattern = fullfile(myFolder, '*.mat');
    matFiles = dir(filePattern);
    orig_tdivList=zeros(1,k);
    change_tdivList=zeros(1,k);
    orig_tchromList=zeros(1,k);
    change_tchromList=zeros(1,k);
    tdivDifList=zeros(1,k);
    tchromDifList=zeros(1,k);
%     unperturbedParams=zeros(1,k);
    counter=0;
    Samples=randsample(length(matFiles),k);
    for i=length(Samples):-1:1
%         try
        tic
        file=Samples(i);
        baseFileName = matFiles(file).name;
        fullFileName = fullfile(myFolder, baseFileName);
        fprintf(1, 'Now reading %s\n', fullFileName);
        matData = load(fullFileName);
        y0=matData.y0;
%         if matData.Sbest >600
%             continue
%         end
        try
            Pset=matData.Pbest;
        catch
            Pset=matData.Pset;
        end
%         Pset(57:59)=[0.24,0.08,3];
%         Pset(8)=0.2;
%         Pset(28)=0;
        [tout, teout, ieout, yout, ~] = CauloSim(celltype,'WT', 1600,Pset,y0.');
        [t,~,cellarrest,ChromRept] = getLastCycleAndStuff(tout,yout,teout,ieout);
        Tdiv1=t(end);
        Tchrom1=ChromRept;
        if strcmp(celltype,'ST')
            if Tchrom1>Tdiv1/2
                Tchrom1=Tchrom1-Tdiv1;
            elseif Tchrom1<=Tdiv1/2
                Tchrom1=Tchrom1;
            end
        end
        if cellarrest == true
            disp('original cell cycle arrested')
            orig_tdivList(i)=[];
            change_tdivList(i)=[];
            orig_tchromList(i)=[];
            change_tchromList(i)=[];
            tdivDifList(i)=[];
            tchromDifList(i)=[];
            continue
        end
        try
            Pset(paramChangesStand(:,1))=paramChangesStand(:,2);
        catch
        end
        try
            Pset(paramChangesMult(:,1))=paramChangesMult(:,2).'.*Pset(paramChangesMult(:,1));
        catch
        end
        [tout, teout, ieout, yout, ~] = CauloSim(celltype,'WT', 1600,Pset,y0.');
        [t,y,cellarrest,ChromRept] = getLastCycleAndStuff(tout,yout,teout,ieout);
        Tdiv2=t(end);
        Tchrom2=ChromRept;
        if strcmp(celltype,'ST')
            if Tchrom2>Tdiv2/2
                Tchrom2=Tchrom2-Tdiv2;
            elseif Tchrom2<=Tdiv2/2
                Tchrom2=Tchrom2;
            end
        end
        if cellarrest == true
            disp('pertrubed cell cycle arrested')
            orig_tdivList(i)=[];
            change_tdivList(i)=[];
            orig_tchromList(i)=[];
            change_tchromList(i)=[];
            tdivDifList(i)=[];
            tchromDifList(i)=[];
            counter=counter+1;
            continue
        end
        
        orig_tdivList(i)=Tdiv1;
        change_tdivList(i)=Tdiv2;
        orig_tchromList(i)
        Tchrom1
        orig_tchromList(i)=Tchrom1;
        change_tchromList(i)
        Tchrom2
        change_tchromList(i)=Tchrom2;
        tdivDifList(i)=Tdiv2-Tdiv1;
        tchromDifList(i)=Tchrom2-Tchrom1;
        toc
%         catch
%             orig_tdivList(i)=[];
%             change_tdivList(i)=[];
%             orig_tchromList(i)=[];
%             change_tchromList(i)=[];
%             tdivDifList(i)=[];
%             tchromDifList(i)=[];
%             counter=counter+1;
%         end
    end
    %save in new 'EvaluatingParams' Folder
    %save std and mean of shifts in tdiv and tchrom. 
    %save arrestcounter
    %save total paramsets evaluated
    %save each param perturbation in a different file.
    %saved file name should have the strain, parameter=?, and cell type
%     try
        nameofFile=strcat(Name,'_',CtrAExpProf,'_',celltype);
        nameofFile=strcat(nameofFile,'.mat');
%     catch
%         nameofFile='tempname.mat';
%     end
    fullFileName = fullfile('ParamCatalog/ParamChangeSims', nameofFile);
    TotalEvals=k-counter;
    ArrestFraction=counter/(length(tdivDifList)+counter);
    WTarrestFraction=(k-(length(tdivDifList)+counter))/k;
    save(fullFileName,'TotalEvals','change_tchromList','orig_tchromList','change_tdivList','orig_tdivList','tchromDifList','tdivDifList','ArrestFraction','WTarrestFraction')
end

