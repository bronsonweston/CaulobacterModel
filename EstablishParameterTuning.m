close all
clc
% clear all
global checkmeouty0 
type='Quick';
CtrAversion='Complete';
celltype = 'SW';
% load('G:\My Drive\Caulobacter Project\CauloSTSW\Parameter Estimation\ParamCatalog\CtrABindingParams_All\01-Jul-2020 13.33.29 Svalue_13494.3688.mat')
% fname='ParamCatalog/SlowParams_All/28-Jun-2020 03.06.02 Svalue_759.2577.mat';
% fname='10-Oct-2020 11.36.49 Svalue_550.861.mat';
% load('ParamCatalog\SlowParams_All\28-Jun-2020 03.06.02 Svalue_759.2577.mat')
% load('ParamCatalog\CompleteParams_All\01-Jul-2020 00.02.42 Svalue_552.9536.mat')
% w = warning('query','last');
% id= w.identifier;F
% warning('off',id);
% warning('on','all');
% seedfile = findSeedParamSet(type,fname);
% return
% parameterstuned = [1:33,35:105, 107:120, 122:130,136];
%parameterstuned = [1:33,35:39, 41:65, 67:82, 85:105, 107:120, 122:130,132,136,142];
parameterstuned = [1:10,12:33,35:39, 41:65, 67:82, 85:105, 107:120, 122:130,136];
% parameterstuned = [1:3,6,8:10, 12:33,35:39, 41:45,47,49:65, 67:82, 85:105, 107:120, 122:130,132,136,142];
load('paramDictionary.mat','paramDictionary');
%%
Pset=Pbest;


%% EvalParamSet
evalparamset=true;
% Pset(3)=0.005;
% Pset(131)=0;
[tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype,'WT', 1600,Pset,y0.');
y0=checkmeouty0;
% mutant='PpleC::Tn & divL(Y550F)';
% [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype,mutant, 1600,Pset,y0.');

if evalparamset == true
    tic
    [Failures] = EvaluateMutants(Pset,y0,1)
    toc
end


% mutant='cckA(Y514D) & PdivK::Tn';
%mutant='deltapopA & deltaPleD';
% mutant='deltaGcrA';
% % % mutant='LS1';
% mutant='SM921';
% Pset(46)=Pset(46)/4;
% Pset(61)=Pset(61)/1.25;
% Pset(89)=Pset(89)/3;
% mutant='PpleC::Tn & deltadivJ';
% mutant='DNAdamage';
% mutant='cckA(Y514D) & PdivK::Tn';
% mutant='LS1';
% mutant='PpleC::Tn & divL(Y550F)';
% mutant='EtOH_Stress';
% mutant='SM921';
% Pset(131)=0;
% Pset(46)=Pset(46)/4;
% [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype,mutant, 1600,Pset,y0.');
% graphCellCycle(tout, teout, ieout, yout, cellCycleIniTimes, Pset, celltype)
% Pset(139)=0;
% [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype,mutant, 1600,Pset,y0.');

% FindCompleteCellCyclesScore({'WT','PpleC::Tn','PdivK::Tn','divKxyl','deltaPleD', 'deltaccrM','RiboActrA-DD', 'deltaPdeA','CtrAD51Edelta3omega','SM921','deltaGcrA','cdG0','divL(A601L)','divL(Y550F)','cckA(Y514D)','deltapopA','deltadgcB'},Pset)

% mutantlist={'WT','PpleC::Tn','PdivK::Tn','deltaPleD', 'deltaccrM','CtrAdelta3omega', 'deltaPdeA','divKxyl','deltaGcrA','deltapopA','divL(A601L)','DNAdamage'};
% mutantlist={'WT','PpleC::Tn','PdivK::Tn','deltaPleD', 'deltaccrM','RiboActrA-DD', 'deltaPdeA','divKxyl','deltaGcrA','deltapopA','divL(A601L)','DNAdamage'};
% mutantlist={'WT','PpleC::Tn', 'deltaccrM', 'deltaPdeA','deltaGcrA','deltapopA','deltapopA & deltaPleD','PdivK::Tn','divL(A601L)','LS1','cckA(Y514D) & PdivK::Tn','deltadgcB'};
% mutantlist={'WT','PpleC::Tn', 'CtrAdelta3omega', 'deltaccrM', 'deltaPdeA','deltaGcrA','deltapopA & PdivK::Tn','deltapopA & deltaPleD','divL(A601L)','LS1'};
% mutantlist={'WT','PpleC::Tn', 'CtrAdelta3omega', 'deltaccrM','deltaGcrA','deltapopA & deltaPleD','divL(A601L)','LS1'};
% mutantlist={'WT','PpleC::Tn', 'CtrAdelta3omega', 'deltaccrM','deltaGcrA','deltapopA & PdivK::Tn','divL(A601L)','LS1'};
% mutantlist={'WT','PpleC::Tn', 'CtrAdelta3omega', 'deltaccrM','deltaGcrA','deltapopA & PdivK::Tn','divL(A601L)','LS1','deltaPleD','deltaSpmX'};
% mutantlist={'WT','deltaPdeA','PpleC::Tn', 'deltapopA & PdivK::Tn','deltaccrM','deltaGcrA','deltapopA & deltaPleD','divL(A601L)','LS1','deltaSpmX','CtrAdelta3omega'};
mutantlist={'WT','deltaPdeA','PpleC::Tn', 'deltapopA & PdivK::Tn','deltaccrM','deltaGcrA','deltaPleD','divL(A601L)','LS1','CtrAdelta3omega','cdG0','CtrAD51E'};
%cckA(Y514D),'deltapopA & deltaPleD'
% mutantlist={'WT','deltaGcrA'};


%% Param Estimation
paramEstimation= false;
% CtrAversion='CtrAbindingOff';
% CtrAversion='Complete';
% CtrAversion='Slow';

% CtrAversion='Partial';
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
    for i=1:10
        [Pset1,~,~]=parameterEstimationMCMCsamples(10,Pset,parameterstuned,5,6,600,1,dirname, mutantlist, stalkon,CtrAversion);
        [Pset,~,~]=parameterEstimationMCMCsamples(5,Pset1,parameterstuned,10,6,250,1,dirname, mutantlist, stalkon,CtrAversion);
%         [Pset,~,~]=parameterEstimationMCMCsamples(10,Pset1,parameterstuned,10,10,100,1,dirname, mutantlist, stalkon,CtrAversion);
    end
end


% cdGT=yout(:, 31)+yout(:, 18)+2*(yout(:, 38)+yout(:, 39)+yout(:, 41)+yout(:, 43));
% cdGavg=trapz(tout,cdGT)/tout(end)

% mutantlist={'WT','divL(A601L)','divL(Y550F)','RiboActrA-DD'};
% Pset(134)=-1*log(2/0.92)/cellCycleIniTimes(end);
% [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype,'LS1', 1500,Pset,y0.');
% [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim('SW','CtrAdelta3omega', 1600,Pset,y0.');


% try
%    graphCellCycle(tout, teout, ieout, yout, cellCycleIniTimes, Pset, celltype)
   graphCellCycle2(tout, teout, ieout, yout, cellCycleIniTimes, Pset, celltype,type)
%    [t,y,cellarrest,ChromRept,Zt] = getLastCycleAndStuff(tout,yout,teout,ieout);
%    ChromRept
%    Zt
%     graphCellCycle3(tout, teout, ieout, yout, cellCycleIniTimes, Pset, celltype)
%     [t,y,cellarrest] = getLastCycle(tout,yout,teout,ieout);
%     figure()
%     plot(t,y(:,11)./(y(:,11)+y(:,12)));
%     plot(t,y(:,1)./(y(:,1)+y(:,2)));
%     pks=findpeaks(y(:, 39));
%     pksT = find(y(:,39)==pks);
%     methCtrA=Pset(136); h_ctra= y(:,15); ksCtrAP1=Pset(1); JiCtrACtrA = Pset(5);
%     CtrA = y(:,1); CtrAP=y(:,2); thetaCtrA = Pset(130); GcrA = y(:,4);
%     JaCtrAGcrA= Pset(6); ksCtrAP2=Pset(2) ; JaCtrACtrA= Pset(7);
%     CtrAP1 = (1-methCtrA+methCtrA.*(2.*(1-h_ctra))).*ksCtrAP1.*(JiCtrACtrA^2./(JiCtrACtrA^2 + (CtrA+CtrAP).^2)).*(1-thetaCtrA+thetaCtrA*GcrA./(JaCtrAGcrA + GcrA));
% 
%     CtrAP2 = ksCtrAP2*(CtrA+CtrAP).^2./(JaCtrACtrA^2 + (CtrA+CtrAP).^2);
%     figure()
%     plot(t,CtrAP1)
%     title('CtrAP1')
%     figure()
%     plot(t,CtrAP2)
% %     title('CtrAP2')
% catch
% end


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

%%
% close all
% [t, y, yayornay, events]=getLastCycle2(tout,yout,teout,ieout);
% ConcentrationRings(t,y,[2,4,48,8],{'CtrA~P','GcrA','DnaA','CcrM'},[1 0 0; 0 0.9 0; 0 0 0.9; 0.9 0 0.9], 3.5, 3.75,events)
% ConcentrationRings(t,y,[2,48,4,8],{'CtrA~P','DnaA','GcrA','CcrM'},[1 0 0; 0 0 0.9; 0 0.9 0; 0.9 0 0.9], 3.5, 3.75,events)

% Pset=paramA;
% parameterstuned = 'All';
% for i=1:10
%     param=Pset;
%     [Pset,BestStorage,ScoreStorage]=parameterEstimationMCMC(param,parameterstuned,5,100,50,false, {'PpleC::Tn'});
%     close all
%     figure()
%     plot(1:length(ScoreStorage),ScoreStorage);
%     hold on
%     scatter(BestStorage,ScoreStorage(BestStorage),'x')
%     addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
% %     ylim([0 2])    
%     [Pset,BestStorage,ScoreStorage]=parameterEstimationMCMC(Pset,parameterstuned,2.5,100,100,true, {'PpleC::Tn'});
%     [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim('SW','PpleC::Tn', 600,Pset);
%     graphCellCycle(tout, teout, ieout, yout, cellCycleIniTimes, param, 'SW')
%     figure()
%     plot(1:length(ScoreStorage),ScoreStorage);
%     hold on
%     scatter(BestStorage,ScoreStorage(BestStorage),'x')
%     addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
% %     ylim([0 2])
%     
% end
% 
% [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim('SW','WT', 800,Pset);
% graphCellCycle(tout, teout, ieout, yout, cellCycleIniTimes, param, 'SW')

%%***********************************************************************%%



% [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype,'WT', 1800,Pset);
% graphCellCycle(tout, teout, ieout, yout, cellCycleIniTimes, Pset, 'ST')
% [t, y, yayornay]=getLastCycle(tout,yout,teout,ieout);
% ConcentrationBars(t,y,[2,3,4,8],{'CtrA~P','GcrA','DnaA','CcrM'},[1 0 0; 0 1 0; 0 0 1; 1 0 1])
% 
% try
%     ScoreStorage=ScoreStorage;
%     figure()
%     plot(1:length(ScoreStorage),ScoreStorage);
%     hold on
%     scatter(BestStorage,ScoreStorage(BestStorage),'x')
%     addToolbarExplorationButtons(gcf) % Adds buttons to figure toolbar
%     ylim([0 400])
% catch
% end
% 

% for i=1:length(yout(1,:))
%     if any(yout(:,i)<=0)
%         i
%     end
% end
