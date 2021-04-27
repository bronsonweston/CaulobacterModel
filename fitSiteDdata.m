% This script plots our Cori site d occupancy fit to data from Siam and Marczynski 2000
% Weston et al. 2021, Cell Systems
close all
clear all 
clc
%CtrA (unphosphorylated) binding with mut-d
S1 = [0 0.25 0.56 1.25 2.25]; 
I1 = [1 0.775 0.55 0.46 0.39];
kd1=2*0.5324;
kd2=0.07752/2;
kd3=0.00017/2;
hFig = figure();
set(hFig, 'Position', [100 100 525 600])

subplot(3,1,1)
scatter(S1,I1,'rs')
hold on
S=0:0.1:2.5;
plot(S, kd1./(kd1+S), 'k','linewidth',1.1)
ylabel('Intensity')
xlabel('[CtrA_U] (\muM)')
title('A')
h = get(gca, 'Title');
set(h, 'Units', 'normalized');
set(h, 'Position', [-0.10,1.15,0]);
set(gca,'FontSize', 12,'fontname', 'Times New Roman');
box on

%CtrA (unphosphorylated) binding no mutation

S2= [0 0.225 0.56 1.25 2.25];
I2= [1 0.35 0.25 0.225 0.15];
I=@(x) ((kd1+x)./(kd1+2*x+(x.^2./kd2)));
subplot(3,1,2)
scatter(S2,I2,'rs')
hold on
S=0:0.05:2.5;
plot(S, I(S), 'k','linewidth',1.1)
ylabel('Intensity')
xlabel('[CtrA_U] (\muM)')
title('B')
ylim([-0.05 1])
h = get(gca, 'Title');
set(h, 'Units', 'normalized');
box on
set(h, 'Position', [-0.10,1.15,0]);
set(gca,'FontSize', 12,'fontname', 'Times New Roman');

%CtrA~P (phosphorylated) binding no mutation
S3=0.006./I2-0.006;
subplot(3,1,3)
scatter(S3, I2, 'rs')
hold on
I=@(x) ((kd1+x)./(kd1+2*x+(x.^2./kd3)));
dots=0:0.001:0.035;
plot(dots, I(dots), 'k','linewidth',1.1)
ylabel('Intensity')
xlabel('[CtrA~P] (\muM)')
title('C')
h = get(gca, 'Title');
set(h, 'Units', 'normalized');
set(h, 'Position', [-0.10,1.15,0]);
set(gca,'FontSize', 12,'fontname', 'Times New Roman');
box on
print(hFig,['Figures/','CtrAbindingFit.png'],'-dpng','-r300');  