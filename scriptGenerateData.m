%% Generate All Maps
%  Maps Should be in the following order:
%                   SlowMap: [72×1 containers.Map]
%                   CompleteMap: [72×1 containers.Map]
%                   ctrAcoriMap: [72×1 containers.Map]
%                   Slow_CoriMut_Map: [72×1 containers.Map]
%                   ctrACori_restoredbindingMap: [72×1 containers.Map]
%                   Complete_CoriMut_Map: [72×1 containers.Map]
%                   Slow_MystKin: [72×1 containers.Map]
%                   Quick_MystKin: [72×1 containers.Map]


evalListofParams(150)
establishNewMap('Slow_CoriMut_Map','SlowParams_All',[131, 0],150);
establishNewMap('ctrACori_restoredbindingMap','CtrABindingParams_All',[131, 1],150);
establishNewMap('Quick_CoriMut_Map','CompleteParams_All',[131, 0],150);
establishNewMap('Slow_MystKin','SlowParams_All',[16, 0.005],150);
establishNewMap('Quick_MystKin','CompleteParams_All',[16, 0.005],150);

%% Get Data on G1/G2 Arrest
% cellcycleViabilityTable('SW','SlowParams_All',150,[],[],'Slow_Normal_SW')
% cellcycleViabilityTable('ST','SlowParams_All',150,[],[],'Slow_Normal_ST')
% cellcycleViabilityTable('SW','SlowParams_All',150,[131, 0],[],'Slow_Corioff_SW')
% cellcycleViabilityTable('ST','SlowParams_All',150,[131, 0],[],'Slow_Corioff_ST')
% 
% cellcycleViabilityTable('SW','CompleteParams_All',150,[],[],'Quick_Normal_SW')
% cellcycleViabilityTable('ST','CompleteParams_All',150,[],[],'Quick_Normal_ST')
% cellcycleViabilityTable('SW','CompleteParams_All',150,[131, 0],[],'Quick_Corioff_SW')
% cellcycleViabilityTable('ST','CompleteParams_All',150,[131, 0],[],'Quick_Corioff_ST')


%% Generate all "ChangingBindingSims"

parfor i=1:6
    if i==1
        removeCoriBinding('SW','SlowParams_All',150)
    elseif i== 2
        removeCoriBinding('ST','SlowParams_All',150)
    elseif i== 3
        removeCoriBinding('SW','CompleteParams_All',150)
    elseif i== 4
        removeCoriBinding('ST','CompleteParams_All',150)
    elseif i== 5
        removeCoriBinding('SW','CtrABindingParams_All',150)
    elseif i== 6
        removeCoriBinding('ST','CtrABindingParams_All',150)
    end
end

%% Generate All "ParamChangeSims"

% mutantParamEval('ST','CompleteParams_All',100,[],[101,1e-4;22,0.1;23 0.1],'deltapopA_deltapled_PdivKTn')
parfor i=1:4
    if i==1
        mutantParamEval('SW','SlowParams_All',150,[46,0;61,0],[22,0.25;23, 0.25],'deltapopA_deltapled_PdivKTn')
    elseif i== 2
        mutantParamEval('ST','SlowParams_All',150,[46,0;61,0],[22,0.25;23, 0.25],'deltapopA_deltapled_PdivKTn')
    elseif i==3
        mutantParamEval('SW','CompleteParams_All',150,[46,0;61,0],[22,0.25;23, 0.25],'deltapopA_deltapled_PdivKTn')
    elseif i==4
        mutantParamEval('ST','CompleteParams_All',150,[46,0;61,0],[22,0.25;23, 0.25],'deltapopA_deltapled_PdivKTn')
    end
end

parfor i=1:4
    if i==1
        mutantParamEval('SW','SlowParams_All',150,[],[101,1e-4;22,0.1;23 0.1],'cckA(Y514D)_PdivKTn')
    elseif i== 2
        mutantParamEval('ST','SlowParams_All',150,[],[101,1e-4;22,0.1;23 0.1],'cckA(Y514D)_PdivKTn')
    elseif i==3
        mutantParamEval('SW','CompleteParams_All',150,[],[101,1e-4;22,0.1;23 0.1],'cckA(Y514D)_PdivKTn')
    elseif i==4
        mutantParamEval('ST','CompleteParams_All',150,[],[101,1e-4;22,0.1;23 0.1],'cckA(Y514D)_PdivKTn')
    end
end

parfor i=1:10
    if i==1
        mutantParamEval('SW','CompleteParams_All',150,[],[4 1.5],'clpxp1.5')
    elseif i==2
        mutantParamEval('SW','CompleteParams_All',150,[131 0],[4,1.5],'clpxp1.5andcor1off')
    elseif i==3
        mutantParamEval('SW','CompleteParams_All',150,[],[4 1.25],'clpxp1.25')
    elseif i==4
        mutantParamEval('SW','CompleteParams_All',150,[131 0],[4,1.25],'clpxp1.25andcor1off')
    elseif i==5
        mutantParamEval('SW','CompleteParams_All',150,[],[4 0.75],'clpxp0.75')
    elseif i==6
        mutantParamEval('SW','CompleteParams_All',150,[131 0],[4,0.75],'clpxp0.75andcor1off')
    elseif i==7
        mutantParamEval('SW','CompleteParams_All',150,[],[4 0.5],'clpxp0.5')
    elseif i==8
        mutantParamEval('SW','CompleteParams_All',150,[131 0],[4,0.5],'clpxp0.5andcor1off')
    elseif i==9
        mutantParamEval('SW','CompleteParams_All',150,[],[4 2],'clpxp2')
    elseif i==10
        mutantParamEval('SW','CompleteParams_All',150,[131 0],[4,2],'clpxp2andcor1off')
    end
end
%%
getCtrAlevels('SW')