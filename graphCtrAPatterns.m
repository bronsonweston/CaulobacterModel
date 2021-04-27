function [] = graphCtrAPatterns()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Holtz=[0.0395516091079948,65.5422682159037;19.7635917140056,9.06182129105349;40.2805696712215,0.816500572097269;59.6832734452962,1.64103264399287;79.3772605120592,45.3412324544610;99.8572501146480,91.1027624446681;119.509625282456,94.8131567681984;139.617419566609,64.7177361440081];
Holtz(:,2)=Holtz(:,2)/100;
Chen = [0 1; 20	0.443745329; 40	0.030918983 ;60	0.009051141;80	0.039503093;100	0.474860095;120	0.695567105; 140	0.70876683];
Mcgrath = [0	12915.912; 20	9480.648; 40	2687.648; 60	515.698; 80	5769.719; 100	10582.426; 120	13602.205;140	13106.439];
Mcgrath(:,2)=Mcgrath(:,2)/max(Mcgrath(:,2));
Domain = [0 1; 20 0.82242751; 40 0.101187217; 60 0.052018273; 80 0.300359313; 100 0.930109248; 120 0.899746897; 140	0.716936547];
Collier = [0 0.735992231; 20 0; 40 0; 60 0; 80 0.457266962; 100 1; 120 0.810964379; 140 0.568760677];
set(0,'DefaultAxesFontName', 'Arial')

figure()
set(gcf,'position',[450 340 341*1.5 142*1.5])
% x = [15 30 30 15];
% y = [0 0 1.2 1.2];
% patch(x,y,'y','FaceAlpha',.5,'EdgeColor','none')
box on
hold on
x = Holtz(:,1);
y = Holtz(:,2);

plot(x,y,'ks','markerfacecolor','k')
hold on
plot(x,y,'k-')
x = Mcgrath(:,1);
y = Mcgrath(:,2);
plot(x,y,'ko','markerfacecolor','k')
plot(x,y,'k-')
x = Domain(:,1);
y = Domain(:,2);
plot(x,y,'k+','markerfacecolor','k','LineWidth',2)
plot(x,y,'k-')
x = Collier(:,1);
y = Collier(:,2);
plot(x,y,'k*','markerfacecolor','k')
plot(x,y,'k-')
ylim([0 1.05])
xlim([0 140])
xlabel('Time (min)')
ylabel('Normalized [CtrA]_T')
set(gca,'FontSize',14)
% print(fig,['Figures/','CtrApatterns.png'],'-dpng','-r300');  

% figure()
% x = Collier(:,1);
% y = Collier(:,2);
% xx = 0:.01:140;
% yy = spline(x,y,xx);
% yy(yy<0)=0;
% plot(x,y,'rs',xx,yy,'r-')
% hold on
% 
% x = Mcgrath(:,1);
% y = Mcgrath(:,2);
% xx = 0:.01:140;
% yy = spline(x,y,xx);
% yy(yy<0)=0;
% plot(x,y,'bs',xx,yy,'b-')
% 
% 
% ylim([0 1.05])

figz=figure()
set(gcf,'position',[450 340 341*1.6 142*1.5])
% x = [15 30 30 15];
% y = [0 0 1.2 1.2];
% patch(x,y,'y','FaceAlpha',.5,'EdgeColor','none')
box on
hold on
plot([18 18],[0 1.1], 'r','linewidth',1.5)
hold on
plot([70 70],[0 1.1], 'r','linewidth',1.5)

x = Mcgrath(:,1);
y = Mcgrath(:,2);
plot(x,y,'ko','markerfacecolor','k')
hold on
plot(x,y,'k-')
x = Collier(:,1);
y = Collier(:,2);
plot(x,y,'kx','markerfacecolor','k', 'markersize',8,'linewidth', 2)
plot(x,y,'k-')
ylim([0 1.05])
xlim([0 140])
xlabel('Time (min)')
ylabel('[CtrA]_T')
set(gca,'FontSize',14)
print(figz,['Figures/','CtrApatterns.png'],'-dpng','-r300');  

fig=figure()
set(gcf,'position',[450 340 341*1.6 142*1.5])
% x = [15 30 30 15];
% y = [0 0 1.2 1.2];
% patch(x,y,'y','FaceAlpha',.5,'EdgeColor','none')
box on
hold on
plot([18 18],[0 1.1], 'r','linewidth',1.5)
hold on
plot([70 70],[0 1.1], 'r','linewidth',1.5)

x = Mcgrath(:,1);
y = Mcgrath(:,2);
plot(x,y,'ko','markerfacecolor','k')
hold on
plot(x,y,'k-')
x = Collier(:,1);
y = Collier(:,2);
plot(x,y,'kx','markerfacecolor','k', 'markersize',8,'linewidth', 2)
plot(x,y,'k-')
x = Holtz(:,1);
y = Holtz(:,2);
plot(x,y,'ks','markerfacecolor','k', 'markersize',8)
plot(x,y,'k-')
x = Domain(:,1);
y = Domain(:,2);
plot(x,y,'k+','markerfacecolor','k','LineWidth',2)
plot(x,y,'k-')
ylim([0 1.05])
xlim([0 140])
xlabel('Time (min)')
ylabel('[CtrA]_T')
set(gca,'FontSize',14)
print(fig,['Figures/','CtrApatterns2.png'],'-dpng','-r300');  
end

