function [] = boxplotScores()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Folders={'ParamCatalog/SlowParams_All','ParamCatalog/CompleteParams_All','ParamCatalog/CtrABindingParams_All'};
l=zeros(1,3);
set(0,'DefaultAxesFontName', 'Arial')

for m=1:length(Folders)
    myFolder=Folders{m}; 
    filePattern = fullfile(myFolder, '*.mat');
    matFiles = dir(filePattern);
    l(m)=length(matFiles);
end
AllScores=zeros(sum(l),1);
count=0;
for m=1:length(Folders)
    myFolder=Folders{m}; 
    filePattern = fullfile(myFolder, '*.mat');
    matFiles = dir(filePattern);
    for i=length(matFiles):-1:1
        count=count+1;
        tic
        baseFileName = matFiles(i).name;
        fullFileName = fullfile(myFolder, baseFileName);
        fprintf(1, 'Now reading %s\n', fullFileName);
        matData = load(fullFileName);
        S=matData.Sbest;
        AllScores(count)=S;
        filetracker{count}=myFolder;
    end
end
% figure()
% boxplot(AllScores,filetracker)
% 
% figure()
Slow=AllScores(1:l(1));
Quick=AllScores((l(1)+1):(l(1)+l(2)));
Cori=AllScores((l(1)+l(2)+1):end);
% bar([1,2,3],[mean(Slow),mean(Quick),mean(Cori)])
% figure()
% set(gcf,'position', [   809   269   129   141])
% boxplot(Slow)
% set(gca,'XTickLabel','Slow')
% xlim([0.8 1.2])
% figure()
% set(gcf,'position', [   809   269   129   141])
% boxplot(Quick)
% set(gca,'XTickLabel','Quick')
% xlim([0.8 1.2])
% 
% figure()
% set(gcf,'position', [   809   269   129   141])
% boxplot(Cori)
% set(gca,'XTickLabel','Cori-')
% curYTick = get(gca,'YTick');
% YLabel = cellstr( num2str( curYTick(:), '%.0f') );
% set(gca, 'YTickLabel', YLabel);
% xlim([0.8 1.2])
% 


figure()
set(gcf,'position', [652   300   328    159])
h1=subplot(1,3,1)
boxplot(Quick,'OutlierSize',3)
xlim([0.86 1.14])
ylim([420 630])
set(gca,'YTickLabel')
yticks([450 525 600])
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
% set(gca,'XTickLabel','Quick')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xlabel(['\fontsize{9}Q','\fontsize{8}UICK'])
ax=gca;
ax.YAxis.FontSize = 8;
% ax.XAxis.FontSize= 10;
pos1=get(gca,'position')

h2=subplot(1,3,2)
boxplot(Slow,'OutlierSize',3)
xlim([0.86 1.14])
ylim([450 1550])
yticks([600 1000 1400])
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
% set(gca,'XTickLabel','Slow')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xlabel(['\fontsize{9}S','\fontsize{8}LOW'])
ax=gca;
ax.YAxis.FontSize = 8;
% ax.XAxis.FontSize= 10;
pos2=get(gca,'position')
set(gca,'position',[pos2(1) pos2(2) pos1(3) pos1(4)])

h3=subplot(1,3,3)
boxplot(Cori,'OutlierSize',3)
xlim([0.86 1.14])
ylim([12250 14250])
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xlabel(['\fontsize{9}C','\fontsize{8}ORI','\fontsize{9}^-'])
yticks([12500 13250 14000])
curYTick = get(gca,'YTick');
YLabel = cellstr( num2str( curYTick(:), '%.0f') );
set(gca, 'YTickLabel', YLabel,'fontsize',8);
ax=gca;
ax.YAxis.FontSize = 8;
% ax.XAxis.FontSize= 10;
pos3=get(gca,'position');
xshift=pos3(1)-(pos2(1)+pos2(3));
set(gca,'position',[pos2(1)+pos1(3)+xshift, pos2(2), pos1(3), pos1(4)])


% set(gca,'position',[0.7329, 0.1205, 0.2078, 0.8045])
% linkaxes([h3,h2,h1])
print(gcf,['Figures/','BoxPlotScores.png'],'-dpng','-r600');  

end

