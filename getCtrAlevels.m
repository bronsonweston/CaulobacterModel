function [] = getCtrAlevels(celltype)
numfiles=150;
if isempty(celltype)
    celltype='SW';
end
% folders=[];
folders={'SlowParams_All','CompleteParams_All','CtrABindingParams_All'};
pCtrAT=zeros(3,numfiles);
pCtrAP=zeros(3,numfiles);
aCtrAT=zeros(3,numfiles);
aCtrAP=zeros(3,numfiles);
allCtrAT=zeros(3,numfiles);
allCtrAP=zeros(3,numfiles);
cdGT=zeros(3,numfiles);
for pcol=1:3
    names={};
    myFolder=strcat('ParamCatalog/',folders{pcol}); filePattern = fullfile(myFolder, '*.mat');
    matFiles = dir(filePattern);
    Samples=randsample(length(matFiles),numfiles);
    for f=1:numfiles
        file=Samples(f);
        baseFileName = matFiles(file).name;
        fullFileName = fullfile(myFolder, baseFileName);
        names{f}=fullFileName;
    end
    count=0;
    parfor i=1:numfiles
        try
        tic
        fprintf(1, 'Now reading %s\n', names{i});
        matData = load(names{i});
        y0=matData.y0;
        Pset=matData.Pbest;
        [tout, teout, ieout, yout, ~] = CauloSim(celltype,'WT', 1600,Pset,y0.');
        [t,y,cellarrest,ChromRept,Zt] = getLastCycleAndStuff(tout,yout,teout,ieout);
        pt=1:find(t==Zt)-1;
        at=find(t==Zt):length(t);
        pCtrAT(pcol,i)=max(y(pt,1)+y(pt,2));
        pCtrAP(pcol,i)=max(y(pt,2));
        aCtrAT(pcol,i)=max(y(at,1)+y(at,2));
        aCtrAP(pcol,i)=max(y(at,2));
        allCtrAT(pcol,i)=max(y(:,1)+y(:,2));
        allCtrAP(pcol,i)=max(y(:,2));
        cdGT(pcol,i)=max(y(:, 31)+y(:, 18)+2*(y(:, 38)+y(:, 39)+y(:, 41)+y(:, 43)));
        toc
        catch
            count=count+1;
        end
    end
end
save(['ParamCatalog/CtrAlevels_',celltype,'.mat'],'pCtrAT','pCtrAP','aCtrAT','aCtrAP','allCtrAT','allCtrAT','allCtrAP', 'cdGT')
end

