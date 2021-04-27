function [] = graphCtrALevels()
close all
load('ParamCatalog/CtrAlevels_ST.mat');
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 11)
% CtrATst=[];
% CtrAPst=[];
for i=1:3
    CtrATst(i)=mean(rmoutliers(allCtrAT(i,:)));
    CtrAPst(i)=mean(rmoutliers(allCtrAP(i,:)));
    cdGst(i)=mean(rmoutliers(cdGT(i,:)));
end
load('ParamCatalog/CtrAlevels_SW.mat');
% CtrATst=[];
% CtrAPst=[];
for i=1:3
    CtrATsw(i)=mean(rmoutliers(allCtrAT(i,:)));
    CtrAPsw(i)=mean(rmoutliers(allCtrAP(i,:)));
    cdGsw(i)=mean(rmoutliers(cdGT(i,:)));
end

fig1=figure();
set(gcf,'position',  [ 200   405   515   261])
CtrAT=[CtrATst;CtrATsw].';
bar(1:3,CtrAT)
str={['\fontsize{9}S','\fontsize{8}LOW'],['\fontsize{9}Q','\fontsize{8}UICK'],['\fontsize{9}C','\fontsize{8}ORI','\fontsize{9}^-']};
set(gca, 'XtickLabel',str, 'Xtick',1:numel(str))
ylabel('Max [CtrA]_T')
Leg=legend('ST','SW'); set(Leg,'position',[0.1598 0.7539 0.1359 0.1245]);
xt=[1:3]-0.225;
yt=CtrAT(:,1)+1.5;
ytxt=num2str(CtrAT(:,1),'%.1f');
text(xt,yt,ytxt,'fontsize',8,'fontweight','bold')
xt=[1:3]+0.065;
yt=CtrAT(:,2)+1.5;
ytxt=num2str(CtrAT(:,2),'%.1f');
text(xt,yt,ytxt,'fontsize',8,'fontweight','bold')
print(gcf,['Figures/','ctrAexpressionA.png'],'-dpng','-r300');  

fig2=figure();
set(gcf,'position',  [ 750   405   515   261])
CtrAP=[CtrAPst;CtrAPsw].';
bar(1:3,CtrAP)
set(gca, 'XtickLabel',str, 'Xtick',1:numel(str))
ylabel('Max [CtrA~P]')
xt=[1:3]-0.225;
yt=CtrAP(:,1)+1.5/2;
ytxt=num2str(CtrAP(:,1),'%.1f');
text(xt,yt,ytxt,'fontsize',8,'fontweight','bold')
xt=[1:3]+0.065;
yt=CtrAP(:,2)+1.5/2;
ytxt=num2str(CtrAP(:,2),'%.1f');
text(xt,yt,ytxt,'fontsize',8,'fontweight','bold')
print(gcf,['Figures/','ctrAexpressionB.png'],'-dpng','-r300');  

fig3=figure();
set(gcf,'position',  [ 200   50   515   261])
r=CtrAP./CtrAT;
bar(1:3,r)
set(gca, 'XtickLabel',str, 'Xtick',1:numel(str))
ylabel('Max [CtrA]_T/Max [CtrA~P]')
xt=[1:3]-0.225;
yt=r(:,1)+1.5*0.8/40;
ytxt=num2str(r(:,1),'%.2f');
text(xt,yt,ytxt,'fontsize',8,'fontweight','bold')
xt=[1:3]+0.065;
yt=r(:,2)+1.5*0.8/40;
ytxt=num2str(r(:,2),'%.2f');
text(xt,yt,ytxt,'fontsize',8,'fontweight','bold')
print(gcf,['Figures/','ctrAexpressionC.png'],'-dpng','-r300');  

fig4=figure();
cdG=[cdGst;cdGsw].';
set(gcf,'position',  [ 700   50   515   261])
bar(1:3,cdG)
set(gca, 'XtickLabel',str, 'Xtick',1:numel(str))
ylabel('Max [cdG]')
xt=[1:3]-0.225;
yt=cdG(:,1)+1.5*0.8/40;
ytxt=num2str(cdG(:,1),'%.2f');
text(xt,yt,ytxt,'fontsize',8,'fontweight','bold')
xt=[1:3]+0.065;
yt=cdG(:,2)+1.5*0.8/40;
ytxt=num2str(cdG(:,2),'%.2f');
text(xt,yt,ytxt,'fontsize',8,'fontweight','bold')
ylim([0 0.8])
print(gcf,['Figures/','ctrAexpressionD.png'],'-dpng','-r300');  

fig=figure
set(gcf, 'Position', [260 170 1030 540]);
box on
h(1)=subplot(2,2,1);
set(gca, 'XtickLabel',str, 'Xtick',1:numel(str))
ylabel('Max [CtrA]_T'); 
h(2)=subplot(2,2,2); 
h(3)=subplot(2,2,3);
h(4)=subplot(2,2,4); 
copyobj(allchild(get(fig1,'CurrentAxes')),h(1)); box on; 
copyobj(allchild(get(fig2,'CurrentAxes')),h(2)); box on
copyobj(allchild(get(fig3,'CurrentAxes')),h(3)); box on;
copyobj(allchild(get(fig4,'CurrentAxes')),h(4)); box on;
% ax = findobj(gcf,'type','axes'); %Retrieve the axes to be copied
% % copyobj([Leg,ax],fig);
end

