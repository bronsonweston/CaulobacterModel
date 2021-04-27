close all
clc
warning('off','all');
%% Load An Initial Parameter Set Here

load('ParamCatalog\SlowParams_All\26-Jun-2020 11.08.00 Svalue_1036.9979.mat')

%% SetUp
parameterstuned = [1:10,12:33,35:39, 41:65, 67:82, 85:105, 107:120, 122:130,136]; %Specifies which parameters should be manipulated
load('paramDictionary.mat','paramDictionary');
Pset=Pbest;



%% Param Estimation
paramEstimation= false;
% CtrAversion='CtrAbindingOff';
% CtrAversion='Complete';
CtrAversion='Slow'; %Specifies if parameterizing to Slow, Quick or Cori- parameter set. Must type in 'Slow', 'CtrAbindingOff', or 'Complete'
mutantlist={'WT','deltaPdeA','PpleC::Tn', 'deltapopA & PdivK::Tn','deltaccrM','deltaGcrA','deltaPleD','divL(A601L)','LS1','CtrAdelta3omega','cdG0','CtrAD51E'};

if paramEstimation == true
    Pset=Pbest;
    if strcmp(CtrAversion,'CtrAbindingOff')
        Pset(131)=0;
    else
        Pset(131)=1;
    end
    count=1;
    direxist=1;
    while direxist == 1
        dirname = strcat(CtrAversion, '_params',num2str(count));
        if ~exist(dirname, 'dir')
            mkdir(dirname)
            direxist=0;
        end
        count=count+1;
    end
    stalkon=true;
    for i=1:5
        [Pset1,~,~]=parameterEstimationMCMCsamples(10,Pset,parameterstuned,5,6,400,1,dirname, mutantlist, stalkon,CtrAversion);
        [Pset,~,~]=parameterEstimationMCMCsamples(5,Pset1,parameterstuned,10,6,200,1,dirname, mutantlist, stalkon,CtrAversion);
%         [Pset,~,~]=parameterEstimationMCMCsamples(10,Pset1,parameterstuned,10,10,100,1,dirname, mutantlist, stalkon,CtrAversion);
    end
end

%% Evaluate Parameter Set
evalparamset=true;
if evalparamset == true
    tic
    [Failures] = EvaluateMutants(Pset,y0,1)
    toc
end

try
    [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim('SW','WT', 1600,Pset,y0.');
    graphCellCycle(tout, teout, ieout, yout, cellCycleIniTimes, Pset, celltype)
catch
end


tic
if evalparamset == true
    [score,scoreset]= getCost(Pset,mutantlist,true,y0,CtrAversion);
else
end
toc
try
    for i=1:length(mutantlist)
        disp(strcat(mutantlist(i), ':'))
        disp(strcat('   ','SW:',' ', num2str(scoreset(i,1:4))))
        disp(strcat('   ', 'ST:', ' ', num2str(scoreset(i+length(mutantlist),1:4))))
    end
    disp(strcat('Score =', ' ',num2str(score)))
catch
end
