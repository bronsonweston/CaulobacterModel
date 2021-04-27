function [cost,costStorage] = getCost(Pcurr, mutants,stalkon, y0,CtrAversion)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
WTweight = 4; PleCweight=1; RiboCtrADDWeight=1500; DivJkoWeight=1500; PdivKTnWeight=1500;

cost=0;
if isempty(mutants)
    mutants= {'WT'};
end
siml = length(mutants);
costStorage = zeros(siml+siml*stalkon,4);
iterationsteps = siml*(1+stalkon);
cdGavg= zeros(1,length(mutants)*(1+stalkon));
DivKPavg=zeros(1,length(mutants)*(1+stalkon));
parfor i=1:1:iterationsteps
    if i> length(mutants)
        celltype = 'ST';
        mutant = mutants(i-length(mutants));
        mutant = string(mutant);
    else
        celltype = 'SW';
        mutant = mutants(i);
        mutant = string(mutant);
    end
%     mutant
    try
        [tout, teout, ieout, yout, ~] = CauloSim(celltype,mutant, 1600,Pcurr,y0.');
    catch
        costStorage(i,:)=[inf,inf,inf];
        continue
    end
%     mutant
%     disp(strcat(mutant,' chromosomes = ', num2str(yout(end,49))));
    numcycles= min(sum(ieout==1),sum(ieout==7));
    [t1,y1,t2,y2,t3,y3,cellarrest]=getLast3Cycles(tout,yout,teout, ieout);

    record=zeros(1,4);
    for cycle=1:3
        if cycle==1
            t=t1; y=y1;
        elseif cycle==2
            t=t2; y=y2;
        elseif cycle==3
            t=t3; y=y3;
        end
        
        if cellarrest == false
            if strcmp(celltype,'SW')
                ChromosomeRepTime=min(t(y(:,13)>0));
            else
                ChromosomeRepTime=min(t(intersect(find(y(:,13)>0),find(y(2:end,47)<=0))));
                if isempty(ChromosomeRepTime)
                    ChromosomeRepTime=min(t(y(:,13)>0));
                end
            end
        end
    try 
        cdGT=y(:, 31)+y(:, 18)+2*(y(:, 38)+y(:, 39)+y(:, 41)+y(:, 43));
        cdGavg(i)=trapz(t,cdGT)/t(end);
%         DivKPavg(i)=(trapz(t,y(:, 6))+2*trapz(t,y(:, 12))+trapz(t,y(:, 36))+trapz(t,y(:, 40)))/t(end);
    catch
        cdGavg(i)=0;
    end
    switch mutant
        case 'WT'
            if  strcmp(celltype, 'ST')
                reptarg = t(end)-10;
                repdiff = min(abs(reptarg-ChromosomeRepTime), ChromosomeRepTime+10);
                record(cycle)= ((repdiff)/4)^2 + ((t(end)-115)/6)^2;
            else
                record(cycle) = WTweight*costFunWT(t,y,Pcurr,CtrAversion)+(y(end,49)-1)*20;
            end
        case 'PpleC::Tn'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/5)^2+((ChromosomeRepTime-20)/4)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/5)^2;
                end
            end
        case ''
            
        case 'DivJko'
            if cellarrest == false
                record(cycle)=100000/t(end);
            else
                record(cycle)=5*(numcycles-1);
            end
        case 'CtrAdelta3omega'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=((1-Pcurr(131))*9+1)*abs((t(end)-145)/5)^2+((ChromosomeRepTime-20)/4)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=((1-Pcurr(131))*9+1)*abs((t(end)-115)/5)^2;
                end
            end
        case 'deltapopA'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/10)^2+((ChromosomeRepTime-20)/8)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/10)^2;
                end
            end
        case 'divKxyl'
            if cellarrest == false
                record(cycle)=50000/t(end);
            else
                record(cycle)=5*(numcycles-1);
            end
        case 'PdivK::Tn'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/5)^2+((ChromosomeRepTime-20)/4)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/5)^2;
                end
            end
        case 'deltapopA & deltaPleD'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/10)^2+((ChromosomeRepTime-20)/8)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/10)^2;
                end
            end    

        case 'deltapopA & PdivK::Tn'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/10)^2+((ChromosomeRepTime-20)/8)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/10)^2;
                end
            end    
        case 'deltaGcrA'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-180)/8)^2;
                end
            end
        case 'divL(Y550F)'
            %                 record(cycle)= PleCweight*costFunPleCko(t,y,WTt,WTy,Pcurr, cellarrest);
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
%                     PleCweight*abs((t(end)-145)/8);
                    record(cycle)=PleCweight*abs((t(end)-200)/10)^2;
                end
            end
        case 'deltaPleD'
            
            if strcmp(CtrAversion, 'Complete')
                dpCtrAT= [0 0.8; 20 0; 20 0; 20 0; 40 0; 60 0; 80 0.5; 100 1; 120 .85; 140 0.7]; %Collier 2006 - DnaA couples DNA... (Normalized based on max expression)
                dpCtrAT(:,1)=dpCtrAT(:,1)*t(end-1)/140;
            else
                %Slow and CtrAbindingOff
                Mcgrath = [0	12915.912; 20	9480.648; 40	2687.648; 60	515.698; 80	5769.719; 100	10582.426; 120	13602.205;140	13106.439];
                Mcgrath(:,2)=Mcgrath(:,2)/max(Mcgrath(:,2));
                dpCtrAT = Mcgrath;
                dpCtrAT(:,1)=dpCtrAT(:,1)*t(end-1)/140;
            end
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    CtrATc=500*findSquares(t,y(:,1)+y(:,2),dpCtrAT,25);
                    record(cycle)=abs((t(end)-145)/5)^2+((ChromosomeRepTime-20)/4)^2+CtrATc;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/5)^2;
                end
            end
        case 'divKcs'
            if cellarrest == false
                record(cycle)=50000/t(end);
            else
                record(cycle)=5*(numcycles-1);
            end
        case 'cdG0'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            end
        case 'deltadgcB'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/5)^2+((ChromosomeRepTime-20)/4)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/5)^2;
                end
            end
        case 'cckA(Y514D)'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/5)^2+((ChromosomeRepTime-20)/4)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/5)^2;
                end
            end
        case 'deltaPodJ'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/5)^2+((ChromosomeRepTime-20)/4)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/5)^2;
                end
            end
        case 'deltaccrM'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-180)/15)^2;
                end
            end
        case 'deltaPdeA'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/5)^2+((ChromosomeRepTime-20)/4)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/5)^2;
                end
            end
         
        case 'CtrAD51E'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            else
                if  strcmp(celltype, 'SW')
                    record(cycle)=abs((t(end)-145)/5)^2+((ChromosomeRepTime-20)/4)^2;
                end
                if  strcmp(celltype, 'ST')
                    record(cycle)=abs((t(end)-115)/5)^2;
                end
            end
        case 'LS1'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            elseif strcmp(celltype, 'SW')
                record(cycle)=((t(end)-160)/10)^2;
            end         
        case 'deltaSpmX'
            if cellarrest == true
                record(cycle)=(150/numcycles)^2;
            end   
        case 'divL(A601L)'
            if cellarrest == false
                record(cycle)=50000/t(end);
            else
                record(cycle)=5*(numcycles-1);
            end
        case 'cckA(Y514D) & PdivK::Tn'
            if cellarrest == false
                record(cycle)=50000/t(end);
            else
                record(cycle)=5*(numcycles-1);
            end
        case 'DNAdamage'
            record(cycle)=sqrt(yout(end,49))-sqrt(2);
    end
    end
    costStorage(i,:)=record;
end

cdGWT=(cdGavg(1)+cdGavg(1+length(mutants)))/2;
if isempty(find(strcmp(mutants,'deltaPleD')==1, 1))
   scorecdGpleD = 0;
else
    index=find(strcmp(mutants,'deltaPleD')==1);
    cdGpleD=(cdGavg(index)+cdGavg(index+length(mutants)))/2;
    scorecdGpleD=(15*abs(cdGpleD/cdGWT-(0.09/0.13)))^2;
end
if isempty(find(strcmp(mutants,'deltaPdeA')==1, 1))
   scorecdGpdeA = 0;
else
    index=find(strcmp(mutants,'deltaPdeA')==1);
    cdGpdeA=(cdGavg(index)+cdGavg(index+length(mutants)))/2;
    scorecdGpdeA=(15*abs(cdGWT/cdGpdeA-(0.13/0.16)))^2;
end

% checkmeouty0
% costStorage
scorecdGWT=10*(25*(cdGWT-0.13))^2/2;
scorecdGpdeA=(100*(mean(cdGpdeA)-0.16))^2;
scorecdGpleD=(100*(mean(cdGpleD)-0.09))^2;
% cdGWT
% mean(cdGpdeA)
% mean(cdGpleD)
% DivKPwt=(DivKPavg(1)+DivKPavg(1+length(mutants)))/2
% if isempty(find(strcmp(mutants,'deltaSpmX')==1, 1))
%    scoreDivKPspmx = 0;
% else
%     index=find(strcmp(mutants,'deltaSpmX')==1);
%     DivKPspmx=(DivKPavg(index)+DivKPavg(index+length(mutants)))/2
%     scoreDivKPspmx=(5*abs(DivKPspmx/DivKPwt-0.23)/0.23)^2;
% end
% if isempty(find(strcmp(mutants,'PpleC::Tn')==1, 1))
%    scoreDivKPplec = 0;
% else
%     index=find(strcmp(mutants,'PpleC::Tn')==1);
%     DivKPplec=(DivKPavg(index)+DivKPavg(index+length(mutants)))/2
%     scoreDivKPplec=(15*abs(DivKPplec/DivKPwt-1.7)/1.7)^2;
% end
% 

for i=1:length(costStorage(:,1))
    costStorage(i,4)=std(costStorage(i,1:3));
    cost= sum(costStorage(i,1:3))/3+costStorage(i,4)+cost;
end
% cost=cost+scoreDivKPplec+scoreDivKPspmx+scorecdGpdeA+scorecdGpleD;
cost=cost+scorecdGWT+scorecdGpdeA+scorecdGpleD;
end
