function [] = barChart(yvalues,err,ylimits,ylab,colors,leglabel,size,fname)
% barChart allows necessary customization not available with MATLAB's bar() function.
% yvalues indicates y values of each bar, err indicates error bars for each bar,
% ylimits details the limits of the y axis, ylab indicates the ylabel, 
% colors indicates the color for each bar, leglabel is the labels for the legend,
% size details the size of the figure and fname is the name of the saved png file
% Weston et al. 2021, Cell Systems

if isempty(size)
    size=[560 200];
end
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 11)
% set(0,'DefaultLegendAutoUpdate','off')

width=0.6; gap=0.1; shift=(width+gap)/2; 
index=[1,2,10,3:9,11,12];
% y=rand(1,12);
y=yvalues;
x=zeros(1,12);
for i=1:length(x)/2
    x([2*i-1,2*i])=[2*i-1-shift,2*i-1+shift];
end

% color=[0 0 1; 1 0 0; 0 1 0];
color=colors;
figure()
set(gcf,'position',[403   246   size])
hold on
for i=1:length(y)
    ind=index(i);
    if i<4
        bar_h1= bar(x(ind), y(ind),'FaceColor',color(i,:),'barwidth',width);
    elseif ind == 12
        bar_h1= bar(x(ind), y(ind),'FaceColor',color(3,:),'barwidth',width);
    elseif rem(ind,2)==1
        bar_h1= bar(x(ind), y(ind),'FaceColor',color(1,:),'barwidth',width);
    else
        bar_h1= bar(x(ind), y(ind),'FaceColor',color(2,:),'barwidth',width);
    end
end

U=zeros(length(y),1);
L=zeros(length(y),1);
for i=1:length(y)
    if y(i)> 0
        U(i)=err(i);
    else
        L(i)=err(i);
    end
end
% legend('WT','Cori^{(-)}','CtrA^{(+)}')
if ~isempty(ylimits)
    set(gca,'ylim',ylimits)
end
ylabel(ylab);
legend(leglabel);
xticks(1:2:11)
xticklabels({['\fontsize{11}S','\fontsize{9}LOW','\fontsize{11} SW'],['\fontsize{11}S','\fontsize{9}LOW','\fontsize{11} ST'],['\fontsize{11}Q','\fontsize{9}UICK','\fontsize{11} SW'], ['\fontsize{11}Q','\fontsize{9}UICK','\fontsize{11} ST'], ['\fontsize{11}C','\fontsize{9}ORI','\fontsize{11}^- SW'], ['\fontsize{11}C','\fontsize{9}ORI','\fontsize{11}^- ST']})
set(gca,'fontsize',11)
er=errorbar(x, y,L,U, '.');
er.Color = 'k';
box on
print(gcf,['Figures/',fname,'.png'],'-dpng','-r600');  

end

