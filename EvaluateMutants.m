function [Failures] = EvaluateMutants(Pset,y0,StalkOn)
%UNTITLED2 Summary of this function goes here

% M = containers.Map('KeyType','char','ValueType','any');
Failures = containers.Map('KeyType','char','ValueType','any');

% [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim('SW','WT', 1600,Pset,y0.');
% y0=checkmeouty0;

celltype = 'SW';
mutantList = {'WT','CtrAdelta3omega','CtrAD51E','CtrAD51Edelta3omega','SM921','deltaGcrA','deltaccrM', ...
    'LS1','PpleC::Tn','divLts','divL(A601L)','divL(Y550F)', 'PpleC::Tn & divL(Y550F)', ...
    'deltadivJ','PpleC::Tn & deltadivJ','deltaSpmX','PpleC::Tn & deltaSpmX', 'PdivK::Tn', ...
    'divKcs','divKxyl','cdG0','PdivK::Tn & cdG0','deltapopA','deltaPleD','deltapopA & PdivK::Tn' ...
    'deltapopA & deltaPleD','deltadgcB','deltaPdeA','deltapopA & deltaPleD & PdivK::Tn', ...
    'cckA(Y514D)','cckA(Y514D) & PdivK::Tn','deltaPodJ','deltahdaA','deltaPleD & PdivK::Tn', 'deltaRcdA'};
ArrestList=[0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0];
% mutantList = {'CtrAdelta3omega','CtrAD51E'};
% ArrestList=[0,0];
C = containers.Map(mutantList,ArrestList);
% [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype,mutant, 1000,Pset,y0.');
arrestList=zeros(length(ArrestList),StalkOn+1);

% for i=1:length(mutantList)
%     if StalkOn== 0
%         DataList(i)=SimStorage();
%     else
%         DataList(i,:)=[SimStorage(),SimStorage()];
%     end
% end
parfor i=1:length(mutantList)
    try
        mutant=char(mutantList(i));
        if strcmp(mutant,'divLts')
            time=2000;
        else
            time=1600;
        end
        [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim('SW',mutant, time,Pset,y0.');
        [tout1,yout1,cellarrest1] = getLastCycle(tout,yout,teout,ieout);
        %     SWcell=SimStorage(mutant,'SW',cellarrest1,tout1,yout1,tout1(end), yout1(end,49));
        if StalkOn == 1
            [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim('ST',mutant, time,Pset,y0.');
            [tout2,yout2,cellarrest2] = getLastCycle(tout,yout,teout,ieout);
            %         STcell=SimStorage(mutant,'ST',cellarrest2,tout2,yout2,tout2(end), yout2(end,49));
        else
            cellarrest2= [];
            %         STcell = [];
        end
        arrestList(i,:)= [cellarrest1,cellarrest2];
        %     DataList(i,:)= [SWcell,STcell];
    catch
        disp(mutant)
        arrestList(i,:)= [abs(ArrestList(i)-1),abs(ArrestList(i)-1)];
    end
end

for i=1:length(mutantList)
    %     mutantList(i)
    %     class(mutantList(i))
    %     arrestList(i,:)
    %     M(char(mutantList(i)))=DataList(i,:);
    for j=1:length(arrestList(i,:))
        %         C(char(mutantList(i)))
        if arrestList(i,j)~= C(char(mutantList(i)))
            if j==1
                celltype='_SW';
            else
                celltype='_ST';
            end
            Failures(strcat(char(mutantList(i)),celltype))=1;
        end
    end
end


end

