close all
clc
clear all 

kd1=0.5324*2;
kd2=0.07752/2;
kd3=0.00017/2;
DNAF=@(CtrA,CtrAP) kd1*(kd1+2*CtrA+2*CtrAP+CtrA^2/kd2+CtrAP^2/kd3+2*CtrA*CtrAP/kd2)^-1;
DNA2CtrAP2=@(CtrA,CtrAP) (DNAF(CtrA,CtrAP)*CtrAP^2)/(kd1*kd3);
DNA2CtrA2=@(CtrA,CtrAP) (DNAF(CtrA,CtrAP)*CtrA^2)/(kd1*kd2);
DNACtrA=@(CtrA,CtrAP) 2*(DNAF(CtrA,CtrAP)*CtrA)/(kd1);
DNACtrAP=@(CtrA,CtrAP) (2*DNAF(CtrA,CtrAP)*CtrAP)/(kd1);
DNACtrACtrAP=@(CtrA,CtrAP) 2*(DNAF(CtrA,CtrAP)*CtrAP*CtrA)/(kd1*kd2);

kappad1=0.2; kappad2=0.006;
AltDNAF=@(CtrA,CtrAP) 1/(1+CtrA/kappad1+CtrAP/kappad2);
AltDNACtrAP=@(CtrA,CtrAP) CtrAP*AltDNAF(CtrA,CtrAP)/kappad2;
AltDNACtrA=@(CtrA,CtrAP) CtrA*AltDNAF(CtrA,CtrAP)/kappad1;

% DNA2CtrAP2= (CtrAT/(R+1))^2./(kd3*(R*CtrAT/(R+1))^2/kd2+


%% Graph 1
ctra=30; ctrap=0;
% DNA=DNAF(ctra,ctrap)
% DNA2CtrAP2i=DNA2CtrAP2(ctra,ctrap)
% DNACtrAi=DNACtrA(ctra,ctrap)
% DNA2CtrA2i=DNA2CtrA2(ctra,ctrap)
ctraplist=[0.1,1,10];
ctralist=0:0.1:30;
DNACtrAP2list = zeros(length(ctraplist),length(ctralist));
DNACtrACtrAPlist = zeros(length(ctraplist),length(ctralist));
DNACtrAPlist = zeros(length(ctraplist),length(ctralist));

for i=1:length(ctralist)
    ctra=ctralist(i);
    for j=1:length(ctraplist)
        ctrap=ctraplist(j);
        DNACtrAP2list(j,i)=DNA2CtrAP2(ctra,ctrap);
        DNACtrAPlist(j,i)=DNACtrAP(ctra,ctrap);
        DNACtrACtrAPlist(j,i)=DNACtrACtrAP(ctra,ctrap);
    end
end


hFig = figure();
set(hFig, 'Position', [100 100 950 300])

subplot(1,3,1) 
plot(ctralist,DNACtrAP2list(1,:),'k','linewidth',3);
hold on
plot(ctralist,DNACtrAP2list(2,:),'r','linewidth',3);
plot(ctralist,DNACtrAP2list(3,:),'b','linewidth',3);
ylabel('[DNA:CtrA~P2]')
xlabel('[CtrA] (\muM)')
title('A')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'fontweight','bold','FontSize', 18,'fontname', 'Arial');
box on

subplot(1,3,2) 
plot(ctralist,DNACtrACtrAPlist(1,:),'k','linewidth',3);
hold on
plot(ctralist,DNACtrACtrAPlist(2,:),'r','linewidth',3);
plot(ctralist,DNACtrACtrAPlist(3,:),'b','linewidth',3);
ylabel('[DNA:CtrA~P:CtrA]')
xlabel('[CtrA] (\muM)')
title('B')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'fontweight','bold','FontSize', 18,'fontname', 'Arial');
box on

subplot(1,3,3) 
plot(ctralist,DNACtrAPlist(1,:),'k','linewidth',3);
hold on
plot(ctralist,DNACtrAPlist(2,:),'r','linewidth',3);
plot(ctralist,DNACtrAPlist(3,:),'b','linewidth',3);
ylabel('[DNA:CtrA~P]')
xlabel('[CtrA] (\muM)')
title('C')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'fontweight','bold','FontSize', 18,'fontname', 'Arial');
box on
legend('[CtrA~P] = 0.1','[CtrA~P] = 1','[CtrA~P] = 10')

hold off

%% Graph 2
ctra=30; ctrap=0;
linethickness=1.5;
fontsize=12;
% DNA=DNAF(ctra,ctrap)
% DNA2CtrAP2i=DNA2CtrAP2(ctra,ctrap)
% DNACtrAi=DNACtrA(ctra,ctrap)
% DNA2CtrA2i=DNA2CtrA2(ctra,ctrap)
ctraplist=[0.1,1,10];
ctralist=0:0.1:30;
DNACtrAP2list = zeros(length(ctraplist),length(ctralist));
DNACtrACtrAPlist = zeros(length(ctraplist),length(ctralist));
DNACtrAPlist = zeros(length(ctraplist),length(ctralist));
DNACtrA2list=zeros(length(ctraplist),length(ctralist));
DNAFlist=zeros(length(ctraplist),length(ctralist));
DNACtrAlist=zeros(length(ctraplist),length(ctralist));
        
for i=1:length(ctralist)
    ctra=ctralist(i);
    for j=1:length(ctraplist)
        ctrap=ctraplist(j);
        DNACtrAP2list(j,i)=DNA2CtrAP2(ctra,ctrap);
        DNACtrAPlist(j,i)=DNACtrAP(ctra,ctrap);
        DNACtrACtrAPlist(j,i)=DNACtrACtrAP(ctra,ctrap);
        DNACtrA2list(j,i)=DNA2CtrA2(ctra,ctrap);
        DNAFlist(j,i)=DNAF(ctra,ctrap);
        DNACtrAlist(j,i)=DNACtrA(ctra,ctrap);
    end
end


hFig = figure();
set(hFig, 'Position', [100 100 950*0.85 600*0.9])

subplot(3,2,1) 
plot(ctralist,DNACtrAP2list(1,:),'k','linewidth',linethickness);
hold on
plot(ctralist,DNACtrAP2list(2,:),'r','linewidth',linethickness);
plot(ctralist,DNACtrAP2list(3,:),'b','linewidth',linethickness);
ylabel('[D:2C~P]')
xlabel('[CtrA] (\muM)')
title('A')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'fontweight','bold','FontSize', fontsize,'fontname', 'Arial');
box on

subplot(3,2,2) 
plot(ctralist,DNACtrA2list(1,:),'k','linewidth',linethickness);
hold on
plot(ctralist,DNACtrA2list(2,:),'r','linewidth',linethickness);
plot(ctralist,DNACtrA2list(3,:),'b','linewidth',linethickness);
ylabel('[D:2C]')
xlabel('[CtrA] (\muM)')
title('B')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'fontweight','bold','FontSize', fontsize,'fontname', 'Arial');
box on
% 
subplot(3,2,3) 
plot(ctralist,DNACtrAPlist(1,:),'k','linewidth',linethickness);
hold on
plot(ctralist,DNACtrAPlist(2,:),'r','linewidth',linethickness);
plot(ctralist,DNACtrAPlist(3,:),'b','linewidth',linethickness);
ylabel('[D:C~P]')
xlabel('[CtrA] (\muM)')
title('C')
set(gca,'fontweight','bold','FontSize', fontsize,'fontname', 'Arial');
box on

subplot(3,2,4) 
plot(ctralist,DNACtrAlist(1,:),'k','linewidth',linethickness);
hold on
plot(ctralist,DNACtrAlist(2,:),'r','linewidth',linethickness);
plot(ctralist,DNACtrAlist(3,:),'b','linewidth',linethickness);
ylabel('[D:C]')
xlabel('[CtrA] (\muM)')
title('D')
set(gca,'fontweight','bold','FontSize', fontsize,'fontname', 'Arial');
box on

subplot(3,2,5) 
plot(ctralist,DNAFlist(1,:),'k','linewidth',linethickness);
hold on
plot(ctralist,DNAFlist(2,:),'r','linewidth',linethickness);
plot(ctralist,DNAFlist(3,:),'b','linewidth',linethickness);
ylabel('[D]_F')
xlabel('[CtrA] (\muM)')
title('E')
set(gca,'fontweight','bold','FontSize', fontsize,'fontname', 'Arial');
box on

subplot(3,2,6) 
plot(ctralist,DNACtrACtrAPlist(1,:),'k','linewidth',linethickness);
hold on
plot(ctralist,DNACtrACtrAPlist(2,:),'r','linewidth',linethickness);
plot(ctralist,DNACtrACtrAPlist(3,:),'b','linewidth',linethickness);
ylabel('[D:C:C~P]')
xlabel('[CtrA] (\muM)')
title('F')
set(gca,'fontweight','bold','FontSize', fontsize,'fontname', 'Arial');
box on


legend('[CtrA~P] = 0.1','[CtrA~P] = 1','[CtrA~P] = 10')



hold off

%% Graph 3
ctra_ctrap_ratio_list=[1,10,25,50];
ctratotlist=0:0.0001:3;
DNACtrAP2list = zeros(length(ctra_ctrap_ratio_list),length(ctratotlist));
DNACtrACtrAPlist = zeros(length(ctra_ctrap_ratio_list),length(ctratotlist));
DNACtrA2list = zeros(length(ctra_ctrap_ratio_list),length(ctratotlist));

for i=1:length(ctratotlist)
    ctratot=ctratotlist(i);
    for j=1:length(ctra_ctrap_ratio_list)
        ctra_ctrap_ratio=ctra_ctrap_ratio_list(j);
        ctrap=ctratot/(1+ctra_ctrap_ratio);
        ctra=ctratot-ctrap;
        DNACtrAP2list(j,i)=DNA2CtrAP2(ctra,ctrap);
        DNACtrA2list(j,i)=DNA2CtrA2(ctra,ctrap);
        DNACtrACtrAPlist(j,i)=DNACtrACtrAP(ctra,ctrap);
    end
end


hFig = figure();
set(hFig, 'Position', [100 100 525*1.5 600*1.5])

subplot(1,3,1) 
plot(ctratotlist,DNACtrAP2list(1,:),'m','linewidth',3);
hold on
plot(ctratotlist,DNACtrAP2list(2,:),'r','linewidth',3);
plot(ctratotlist,DNACtrAP2list(3,:),'b','linewidth',3);
plot(ctratotlist,DNACtrAP2list(4,:),'g','linewidth',3);
ylabel('[DNA:CtrA~P2]')
xlabel('[CtrA]_T (\muM)')
title('A')
DNACtrAP2listA=DNACtrAP2list;
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'FontSize', 14,'fontname', 'Arial');
box on

subplot(1,3,2) 
plot(ctratotlist,DNACtrACtrAPlist(1,:),'m','linewidth',3);
hold on
plot(ctratotlist,DNACtrACtrAPlist(2,:),'r','linewidth',3);
plot(ctratotlist,DNACtrACtrAPlist(3,:),'b','linewidth',3);
plot(ctratotlist,DNACtrACtrAPlist(4,:),'g','linewidth',3);
ylabel('[DNA:CtrA~P:CtrA]')
xlabel('[CtrA]_T (\muM)')
title('B')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'FontSize', 14,'fontname', 'Arial');
box on

subplot(1,3,3) 
plot(ctratotlist,DNACtrA2list(1,:),'m','linewidth',3);
hold on
plot(ctratotlist,DNACtrA2list(2,:),'r','linewidth',3);
plot(ctratotlist,DNACtrA2list(3,:),'b','linewidth',3);
plot(ctratotlist,DNACtrA2list(4,:),'g','linewidth',3);
ylabel('[DNA:CtrA2]')
xlabel('[CtrA]_T (\muM)')
title('C')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'FontSize', 14,'fontname', 'Arial');
box on
legend(strcat('[CtrA:CtrA~P] =',num2str(ctra_ctrap_ratio_list(1))),strcat('[CtrA:CtrA~P] =',num2str(ctra_ctrap_ratio_list(2))),strcat('[CtrA:CtrA~P] =',num2str(ctra_ctrap_ratio_list(3))))


%% Graph 4

ctra_ctrap_ratio_list=[1,10,25,50];
ctraplist=0:0.001:1;
DNACtrAP2list = zeros(length(ctra_ctrap_ratio_list),length(ctraplist));
DNAFlist = zeros(length(ctra_ctrap_ratio_list),length(ctraplist));
DNACtrA2list = zeros(length(ctra_ctrap_ratio_list),length(ctraplist));

for i=1:length(ctraplist)
    ctrap=ctraplist(i);
    for j=1:length(ctra_ctrap_ratio_list)
        ctra_ctrap_ratio=ctra_ctrap_ratio_list(j);
        ctra=ctrap*ctra_ctrap_ratio;
        DNACtrAP2list(j,i)=DNA2CtrAP2(ctra,ctrap);
        DNACtrA2list(j,i)=DNA2CtrA2(ctra,ctrap);
        DNAFlist(j,i)=DNAF(ctra,ctrap);
    end
end


hFig = figure();
set(hFig, 'Position', [100 100 1000*1.3*4/3 241*1.1])

subplot(1,4,1) 
plot(ctraplist,DNACtrAP2list(1,:),'m','linewidth',3);
hold on
plot(ctraplist,DNACtrAP2list(2,:),'r','linewidth',3);
plot(ctraplist,DNACtrAP2list(3,:),'b','linewidth',3);
plot(ctraplist,DNACtrAP2list(4,:),'g','linewidth',3);
set(gca,'FontSize', 14,'fontname', 'Arial');

ylabel('[DNA:CtrA~P_2]','fontsize',18)
xlabel('[CtrA~P] (\muM)','fontsize',18)
% title('A')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
box on

subplot(1,4,2) 
% plot(ctraplist,DNACtrA2list(1,:),'m','linewidth',3);
% hold on
% plot(ctraplist,DNACtrA2list(2,:),'r','linewidth',3);
% plot(ctraplist,DNACtrA2list(3,:),'b','linewidth',3);
% plot(ctraplist,DNACtrA2list(4,:),'g','linewidth',3);
% ylabel('[DNA:CtrA_2]')
% xlabel('[CtrA~P] (\muM)')
plot(ctratotlist,DNACtrAP2listA(1,:),'m','linewidth',3);
hold on
plot(ctratotlist,DNACtrAP2listA(2,:),'r','linewidth',3);
plot(ctratotlist,DNACtrAP2listA(3,:),'b','linewidth',3);
plot(ctratotlist,DNACtrAP2listA(4,:),'g','linewidth',3);
set(gca,'FontSize', 14,'fontname', 'Arial');
ylabel('[DNA:CtrA~P_2]','fontsize',18)
xlabel('[CtrA]_T (\muM)','fontsize',18)
xlim([0 1])
% title('B')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
AX=legend({strcat('[CtrA_U]/[CtrA~P] =',num2str(ctra_ctrap_ratio_list(1))),strcat('[CtrA_U]/[CtrA~P] =',num2str(ctra_ctrap_ratio_list(2))),strcat('[CtrA_U]/[CtrA~P] =',num2str(ctra_ctrap_ratio_list(3))),strcat('[CtrA_U]/[CtrA~P] =',num2str(ctra_ctrap_ratio_list(4)))},'FontSize',14);
% set(AX,'position',[0.6610 0.6259 0.1890 0.300])
set(AX,'position',[0.7316 0.4604 0.1280 0.4736])
box on

SatValues=zeros(1,length(DNACtrAP2list(:,1)));
for i=1:length(DNACtrAP2list(:,1))
    SatValues(i)=DNACtrAP2list(i,end);
end
subplot(1,4,3) 
plot(ctratotlist,DNACtrAP2listA(1,:)/SatValues(1),'m','linewidth',3);
hold on
plot(ctratotlist,DNACtrAP2listA(2,:)/SatValues(2),'r','linewidth',3);
plot(ctratotlist,DNACtrAP2listA(3,:)/SatValues(3),'b','linewidth',3);
plot(ctratotlist,DNACtrAP2listA(4,:)/SatValues(4),'g','linewidth',3);
set(gca,'FontSize', 14,'fontname', 'Arial');
ylabel('\delta_{sat}','fontsize',24)
xlabel('[CtrA]_T (\muM)','fontsize',18)

xlim([0 1])

subplot(1,4,4) 
p1=plot(ctraplist,DNAFlist(1,:),'m','linewidth',3);
hold on
p2=plot(ctraplist,DNAFlist(2,:),'r','linewidth',3);
p3=plot(ctraplist,DNAFlist(3,:),'b','linewidth',3);
p4=plot(ctraplist,DNAFlist(4,:),'g','linewidth',3);
ylabel('[DNA]_F')
xlabel('[CtrA~P] (\muM)')
% title('C')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'FontSize', 14,'fontname', 'Arial');
box on
% AX=legend({strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(1))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(2))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(3))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(4)))},'FontSize',12);
% set(AX,'position',[0.7065 0.5491 0.189 0.38])
set(p1, 'visible', 'off');
set(p2, 'visible', 'off');
set(p3, 'visible', 'off');
set(p4, 'visible', 'off');
set(gca, 'visible', 'off');
% AX=legend({strcat('[CtrA]:[CtrA~P]
% =',num2str(ctra_ctrap_ratio_list(1))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(2))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(3))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(4)))},'FontSize',10);0
print(hFig,['Figures/','CtrA2CtrAPratioA_C.png'],'-dpng','-r300');  


%% Graph 5

ctraplist=0:0.001:0.1;
DNACtrAP2list = zeros(1,length(ctraplist));

for i=1:length(ctraplist)
    ctrap=ctraplist(i);
    ctra=0;
    DNACtrAP2list(i)=DNA2CtrAP2(ctra,ctrap);
end


hFig = figure();
% set(hFig, 'Position', [100 100 525 600])
plot(ctraplist,DNACtrAP2list,'k','linewidth',3);

ylabel('[DNA:CtrA~P2]')
xlabel('[CtrA~P] (\muM)')
title('A')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'FontSize', 14,'fontname', 'Arial');
box on

hFig = figure();
% set(hFig, 'Position', [100 100 525 600])
set(hFig, 'Position', [500 100 333*1.1 241*1.1])
plot(ctraplist*6.022e23*1e-15*1e-6,DNACtrAP2list,'k','linewidth',3);
hold on
num=0.0096*6.022e23*1e-15*1e-6; %5.78
plot([0 num], [0.5 0.5], 'k--','linewidth',1)
plot([num num], [0 0.5],'k--','linewidth',1)
xticks([5.78 20 40])
xticklabels({'5.8','20', '40'})
set(gca,'FontSize', 14,'fontname', 'Arial');
ylabel('[DNA:CtrA~P_2]','FontSize', 18)
xlabel('CtrA~P Molecules','FontSize', 18)
% title('D')
xlim([0 40])
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
print(hFig,['Figures/','CtrA2CtrAPratioD.png'],'-dpng','-r300');  

box on
%% Graph 6

ctra_ctrap_ratio_list=2:0.001:800;
DNACtrAP2list = zeros(1,length(ctra_ctrap_ratio_list));
ctraT=15;
ctrA=ctraT*ctra_ctrap_ratio_list./(1+ctra_ctrap_ratio_list);
ctrAP=ctraT./(1+ctra_ctrap_ratio_list);
for j=1:length(ctra_ctrap_ratio_list)
    DNACtrAP2list(j)=DNA2CtrAP2(ctrA(j),ctrAP(j));
end
syms x
find(DNA2CtrAP2(15-x,x)==0.5)

hFig = figure();
set(hFig, 'Position', [500 100 400 400])

plot(ctra_ctrap_ratio_list,DNACtrAP2list,'k','linewidth',3);
ylabel('[DNA:CtrA~P_2]')
xlabel('[CtrA]:[CtrA~P]')
set(gca,'FontSize', 14,'fontname', 'Arial');

hFig = figure();
set(hFig, 'Position', [100 100 400 400])
plot(ctrAP,DNACtrAP2list,'k','linewidth',3);
ylabel('[DNA:CtrA~P_2]')
xlabel('[CtrA~P]')
set(gca,'FontSize', 14,'fontname', 'Arial');
box on

hFig = figure();
set(hFig, 'Position', [100 100 333*1.1 241*1.1])
plot(ctrAP*6.022e23*1e-15*1e-6,DNACtrAP2list,'k','linewidth',3);
hold on
num=0.688*6.022e23*1e-15*1e-6; %414.3
plot([0 num], [0.5 0.5], 'k--','linewidth',1)
plot([num num], [0 0.5],'k--','linewidth',1)
xticks([414.3 1000 2000])
xticklabels({'414.3', '1000', '2000'})
set(gca,'FontSize', 14,'fontname', 'Arial');

ylabel('[DNA:CtrA~P_2]','FontSize', 18)
xlabel('CtrA~P Molecules','FontSize', 18)
xlim([0 2000])
% title('E')
box on
print(hFig,['Figures/','CtrA2CtrAPratioE.png'],'-dpng','-r300');  

%% Alternative Graph 1

ctra=30; ctrap=0;
% DNA=DNAF(ctra,ctrap)
% DNA2CtrAP2i=DNA2CtrAP2(ctra,ctrap)
% DNACtrAi=DNACtrA(ctra,ctrap)
% DNA2CtrA2i=DNA2CtrA2(ctra,ctrap)
ctraplist=[0.1,1,10];
ctralist=0:0.1:30;
DNACtrAPlist = zeros(length(ctraplist),length(ctralist));
DNAFlist = zeros(length(ctraplist),length(ctralist));
DNACtrAlist = zeros(length(ctraplist),length(ctralist));

for i=1:length(ctralist)
    ctra=ctralist(i);
    for j=1:length(ctraplist)
        ctrap=ctraplist(j);
        DNACtrAPlist(j,i)=AltDNACtrAP(ctra,ctrap);
        DNAFlist(j,i)=AltDNAF(ctra,ctrap);
        DNACtrAlist(j,i)=AltDNACtrA(ctra,ctrap);
    end
end


hFig = figure();
set(hFig, 'Position', [100 100 950 300])

subplot(1,3,1) 
plot(ctralist,DNACtrAPlist(1,:),'k','linewidth',3);
hold on
plot(ctralist,DNACtrAPlist(2,:),'r','linewidth',3);
plot(ctralist,DNACtrAPlist(3,:),'b','linewidth',3);
ylabel('[DNA:CtrA~P]')
xlabel('[CtrA] (\muM)')
title('A')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'fontweight','bold','FontSize', 18,'fontname', 'Arial');
box on

subplot(1,3,2) 
plot(ctralist,DNAFlist(1,:),'k','linewidth',3);
hold on
plot(ctralist,DNAFlist(2,:),'r','linewidth',3);
plot(ctralist,DNAFlist(3,:),'b','linewidth',3);
ylabel('[DNAF]')
xlabel('[CtrA] (\muM)')
title('B')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'fontweight','bold','FontSize', 18,'fontname', 'Arial');
box on

subplot(1,3,3) 
plot(ctralist,DNACtrAlist(1,:),'k','linewidth',3);
hold on
plot(ctralist,DNACtrAlist(2,:),'r','linewidth',3);
plot(ctralist,DNACtrAlist(3,:),'b','linewidth',3);
ylabel('[DNA:CtrA]')
xlabel('[CtrA] (\muM)')
title('C')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'fontweight','bold','FontSize', 18,'fontname', 'Arial');
box on
legend('[CtrA~P] = 0.1','[CtrA~P] = 1','[CtrA~P] = 10')

hold off

%% Graph LastThing

ctra_ctrap_ratio_list=[1,10,25,50];
ctraplist=0:0.001:1;
DNACtrAP2list = zeros(length(ctra_ctrap_ratio_list),length(ctraplist));
DNAFlist = zeros(length(ctra_ctrap_ratio_list),length(ctraplist));
DNACtrA2list = zeros(length(ctra_ctrap_ratio_list),length(ctraplist));

for i=1:length(ctraplist)
    ctrap=ctraplist(i);
    for j=1:length(ctra_ctrap_ratio_list)
        ctra_ctrap_ratio=ctra_ctrap_ratio_list(j);
        ctra=ctrap*ctra_ctrap_ratio;
        DNACtrAP2list(j,i)=DNA2CtrAP2(ctra,ctrap);
        DNACtrA2list(j,i)=DNA2CtrA2(ctra,ctrap);
        DNAFlist(j,i)=DNAF(ctra,ctrap);
    end
end


hFig = figure();
set(hFig, 'Position', [   287   161   884   498])

% subplot(1,2,1) 
plot(ctraplist,DNACtrAP2list(1,:),'m','linewidth',3);
hold on
plot(ctraplist,DNACtrAP2list(2,:),'r','linewidth',3);
plot(ctraplist,DNACtrAP2list(3,:),'b','linewidth',3);
plot(ctraplist,DNACtrAP2list(4,:),'g','linewidth',3);
set(gca,'FontSize', 14,'fontname', 'Arial');

ylabel('[DNA:CtrA~P_2]','fontsize',18)
xlabel('[CtrA~P] (\muM)','fontsize',18)
yticks([0 0.5 1])
yticklabels({'0%','50%','100%'})
% title('A')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
box on

% subplot(1,4,2) 
% % plot(ctraplist,DNACtrA2list(1,:),'m','linewidth',3);
% % hold on
% % plot(ctraplist,DNACtrA2list(2,:),'r','linewidth',3);
% % plot(ctraplist,DNACtrA2list(3,:),'b','linewidth',3);
% % plot(ctraplist,DNACtrA2list(4,:),'g','linewidth',3);
% % ylabel('[DNA:CtrA_2]')
% % xlabel('[CtrA~P] (\muM)')
% plot(ctratotlist,DNACtrAP2listA(1,:),'m','linewidth',3);
% hold on
% plot(ctratotlist,DNACtrAP2listA(2,:),'r','linewidth',3);
% plot(ctratotlist,DNACtrAP2listA(3,:),'b','linewidth',3);
% plot(ctratotlist,DNACtrAP2listA(4,:),'g','linewidth',3);
set(gca,'FontSize', 16,'fontname', 'Arial');
ylabel('[DNA:CtrA~P_2]','fontsize',18)
xlabel('[CtrA]_T (\muM)','fontsize',18)
xlim([0 1])
% title('B')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
AX=legend({strcat('[CtrA_U]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(1))),strcat('[CtrA_U]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(2))),strcat('[CtrA_U]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(3))),strcat('[CtrA_U]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(4)))},'FontSize',14,'Location','eastoutside');
% set(AX,'position',[0.6610 0.6259 0.1890 0.300])
% set(AX,'position',[0.7316 0.4604 0.1280 0.4736])
box on
print(hFig,['Figures/','pptBinding.png'],'-dpng','-r300');  

return
% SatValues=zeros(1,length(DNACtrAP2list(:,1)));
% for i=1:length(DNACtrAP2list(:,1))
%     SatValues(i)=DNACtrAP2list(i,end);
% end
% subplot(1,4,3) 
% plot(ctratotlist,DNACtrAP2listA(1,:)/SatValues(1),'m','linewidth',3);
% hold on
% plot(ctratotlist,DNACtrAP2listA(2,:)/SatValues(2),'r','linewidth',3);
% plot(ctratotlist,DNACtrAP2listA(3,:)/SatValues(3),'b','linewidth',3);
% plot(ctratotlist,DNACtrAP2listA(4,:)/SatValues(4),'g','linewidth',3);
% set(gca,'FontSize', 14,'fontname', 'Arial');
% ylabel('\delta_{sat}','fontsize',24)
% xlabel('[CtrA]_T (\muM)','fontsize',18)

xlim([0 1])

subplot(1,4,4) 
p1=plot(ctraplist,DNAFlist(1,:),'m','linewidth',3);
hold on
p2=plot(ctraplist,DNAFlist(2,:),'r','linewidth',3);
p3=plot(ctraplist,DNAFlist(3,:),'b','linewidth',3);
p4=plot(ctraplist,DNAFlist(4,:),'g','linewidth',3);
ylabel('[DNA]_F')
xlabel('[CtrA~P] (\muM)')
% title('C')
% h = get(gca, 'Title');
% set(h, 'Units', 'normalized');
% set(h, 'Position', [-0.10,35,0]);
set(gca,'FontSize', 14,'fontname', 'Arial');
box on
% AX=legend({strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(1))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(2))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(3))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(4)))},'FontSize',12);
% set(AX,'position',[0.7065 0.5491 0.189 0.38])
set(p1, 'visible', 'off');
set(p2, 'visible', 'off');
set(p3, 'visible', 'off');
set(p4, 'visible', 'off');
set(gca, 'visible', 'off');
% AX=legend({strcat('[CtrA]:[CtrA~P]
% =',num2str(ctra_ctrap_ratio_list(1))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(2))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(3))),strcat('[CtrA]:[CtrA~P] =',num2str(ctra_ctrap_ratio_list(4)))},'FontSize',10);0
% print(hFig,['Figures/','CtrA2CtrAPratioA_C.png'],'-dpng','-r300');  
