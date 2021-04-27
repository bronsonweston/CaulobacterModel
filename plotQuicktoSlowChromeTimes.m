function [] = plotQuicktoSlowChromeTimes()
% close all
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 11)
set(0,'DefaultLegendAutoUpdate','off')
kdmultiplier =[0.2,0.4,0.6,0.8,1];
QuickSet= [];
CoriSet= [];
Quickerr=[];
Corierr=[];
for i=1:length(kdmultiplier)
    if kdmultiplier(i)~=1
        file1=['ParamCatalog/ParamChangeSims/clpxp',num2str(kdmultiplier(i)),'_CompleteParams_All_SW.mat'];
        file2=['ParamCatalog/ParamChangeSims/clpxp',num2str(kdmultiplier(i)),'andcor1off_CompleteParams_All_SW.mat'];
    else
        file1='ParamCatalog/ParamChangeSims/clpxp1.0_CompleteParams_All_SW.mat';
        file2='ParamCatalog/ParamChangeSims/clpxp1.0andcor1off_CompleteParams_All_SW.mat';
    end
    load(file1);
    change_tchromList=rmoutliers(change_tchromList);
    QuickSet(i)=mean(change_tchromList);
    Quickerr(i) =std(change_tchromList);
    load(file2);
    change_tchromList=rmoutliers(change_tchromList);
    CoriSet(i)=mean(change_tchromList);
    Corierr(i) =std(change_tchromList);
end

figure()
set(gcf,'position',[   710   246   416   325]);
X=linspace(0,2.5,10);
x=kdmultiplier;
y=QuickSet;
scatter(x,y,'b','filled')
hold on
y=CoriSet;
scatter(x,y,'k','filled')
y=QuickSet;
% m = polyfitn(x,y,1)
[P,S] = polyfit(x,y,1);
yfit = P(1)*X+P(2);
P(1)
disp([num2str(P(1)),'X+',num2str(P(2))])
Rs1= 1 - (S.normr/norm(y - mean(y)))^2
plot(X,yfit,'b--','linewidth',1.5);
errorbar(x,y,Quickerr,'b', 'LineStyle', 'none')

y=CoriSet;
hold on
[P,S] = polyfit(x,y,1);
yfit = P(1)*X+P(2);
plot(X,yfit,'k--','linewidth',1.5);
disp([num2str(P(1)),'X+',num2str(P(2))])
Rs2= 1 - (S.normr/norm(y - mean(y)))^2
errorbar(kdmultiplier,CoriSet,Corierr,'k', 'LineStyle', 'none')


xlim([0.15 1.05])
xlabel('[ClpXP]')
ylabel('t^{cr} (min)')
ylim([15 30])
box on
legend({'WT','WT - CtrA_U:{\it Cori}'})
print(gcf,['Figures/','QuicktoSlowChromTimes.png'],'-dpng','-r600');  
end

