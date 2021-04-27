function [seedfile] = findSeedParamSet(CtrAversioni,fname)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mutantlist={'WT','deltaPdeA','PpleC::Tn', 'deltapopA & PdivK::Tn','deltaccrM','deltaGcrA','deltaPleD','divL(A601L)','LS1','CtrAdelta3omega','cdG0','CtrAD51E'};
parameterstuned = [1:10,12:33,35:39, 41:65, 67:82, 85:105, 107:120, 122:130,136];
CtrAversion=CtrAversioni;
if strcmp('Cori-',CtrAversion)
    CtrAversion=CtrAbindingOff;
elseif strcmp('Quick',CtrAversion)
    CtrAversion='Complete';
end
load(fname,'Pbest','Sbest', 'y0');
% CtrAversion='Partial';
Pset=Pbest;
if strcmp(CtrAversion,'CtrAbindingOff')
    Pset(131)=0;
else
    Pset(131)=1;
end
stalkon=true;
Sbestlast=Pbest;
seedcount=0;
while seedcount<2
    %create directory
    count=1;
    direxist=1;
    while direxist == 1
        dirname = strcat(CtrAversion, '_SeedSearch_params',num2str(count));
        if ~exist(dirname, 'dir')
            mkdir(dirname)
            direxist=0;
        end
        count=count+1;
    end
    %run parameterization algorithm
    [Pset1,~,~]=parameterEstimationMCMCsamples(3,Pset,parameterstuned,35,2.5,250,1,dirname, mutantlist, stalkon,CtrAversion);
    [Pset2,~,~]=parameterEstimationMCMCsamples(20,Pset1,parameterstuned,5,0.7,750,1,dirname, mutantlist, stalkon,CtrAversion);
    [Pset,~,~]=parameterEstimationMCMCsamples(10,Pset2,parameterstuned,15,35,250,1,dirname, mutantlist, stalkon,CtrAversion);
    %find best score from last run
    filePattern = fullfile(dirname, '*.mat');
    matFiles = dir(filePattern);
    for i=length(matFiles):-1:1
        count=count+1;
        tic
        baseFileName = matFiles(i).name;
        fullFileName = fullfile(dirname, baseFileName);
        matData = load(fullFileName);
        S=matData.Sbest;
        if S<Sbest
            Sbest=S;
            seedfile=fullFileName;
        end
    end
    if Sbestlast==Sbest
        seedcount=seedcount+1;
    else
        seedcount=0;
        Sbestlast=Sbest;
    end
end
load(seedfile,'Pbest','Sbest', 'y0');
filename=['seedfile_',CtrAversioni,'.mat'];
save(filename,'Pbest','Sbest', 'y0')
end

