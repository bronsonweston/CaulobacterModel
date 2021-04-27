function [outputArg1,outputArg2] = plotExpressionPatterns(File,k)
% plotExpressionPatterns(File,k) plots wild type, swarmer simulations of
% k parameter sets from the specified File.
% Utilized for Figure 4A&B and Figure 5C of main text
% Weston et al. 2021, Cell Systems

set(0,'DefaultAxesFontName', 'Arial')
gold= [218,165,32]./255;
brown=[165,42,42]./255;
gray=[169,169,169]./255;
close all
dpCpdR = [1*0.93 0.42; 25*0.93 0.938; 50*0.93 0.636; 75*0.93 0.48; 100*0.93 0.07; 125*0.93 0.275; 150*0.93 0]; %Iniesta et al. 2006 Figure 5A estimated
dpcdG = [0 1/3; 20 1; 40 0.55; 60 0.35; 80 0.3; 100 0.25; 120 0.05; 140 0.35]; %Abel 2013
dpcdG(:,2)=dpcdG(:,2)*0.27;
dpSciP = [0 0.87; 20 0.856; 40 0.043; 80 0; 100 0.02; 120 0.4; 140 1]; %Tan et al. 2010
dpCtrAP = [10*1.07 0.4; 80*1.07 0.5; 105*1.07 .8; 130*1.07 1]; %Jacobs 2003 (relative to max = 1) removed second point (t=40, P=0.15) as CtrAT suggests it is incorrect
dpCckAk = [10 0.4; 40 0.3; 70 0.4; 100 1; 130 1]; %Jacobs 2003 (relative to max = 1)
dpGcrA = [0 0; 20 .2; 40 .7;60  1; 80 0.8; 100 0.6; 120 .5; 140 .4];
dpDnaA = [0 0.78;15 1; 30 0.9; 45 0.725; 60 0.375; 75 0.325; 90 0.225; 105 0.42; 120 0.4; 135 0.6]; % Cheng and Keiler 2009 (cell cycle time 136 +- 11 min)
dpDnaA(:,1) = dpDnaA(:,1)*145/136;
dpCcrM = [5 .15; 45 0; 90 0.5; 120 1; 150 0.25]; % Grunenfelder et al. 2001 Proteomic analysis of the bacterial cell cycle (relative to max)
dpCcrM2 = [0 0; 20 0; 40 0; 60 0; 80 0; 100 0.39; 120 1; 140 0.89]; % Zhou et al. 2018; "Cell cycle-controlled clearance of the CcrM..."

if strcmp(File,'CompleteParams_All')
    dpCtrAT = [0 0.8; 20 0; 40 0; 60 0; 80 0.5; 100 1; 120 .85; 140 0.7]; %Collier 2006 - DnaA couples DNA... (Normalized based on max expression)
    dpCtrAT(:,2)=dpCtrAT(:,2)*30;
else
    Mcgrath = [0	12915.912; 20	9480.648; 40	2687.648; 60	515.698; 80	5769.719; 100	10582.426; 120	13602.205;140	13106.439];
    Mcgrath(:,2)=20*Mcgrath(:,2)/max(Mcgrath(:,2));
    dpCtrAT = Mcgrath;
end
CtrATnorm=zeros(1,k);
CtrAPnorm=zeros(1,k);
myFolder=strcat('ParamCatalog/',File); filePattern = fullfile(myFolder, '*.mat');
matFiles = dir(filePattern);
Samples=randsample(length(matFiles),k);
fig=figure(); set(gcf, 'Position', [260 170 453   127]); xlabel('Time (min)'); ylabel('[CtrA]_T (\muM)');
fs=get(gca, 'fontsize'); set(gca, 'fontsize',fs*1.2);

fig1=figure(); xlabel('Time (min)'); ylabel('[CdG] (\muM)');
fig2=figure(); xlabel('Time (min)'); ylabel('[CpdR]');
fig3=figure(); xlabel('Time (min)'); ylabel('[DnaA]');
fig4=figure(); xlabel('Time (min)'); ylabel('[GcrA]');
fig5=figure(); xlabel('Time (min)'); ylabel('[CcrM]');
fig6=figure(); xlabel('Time (min)'); ylabel('[CtrA~P] (\muM)');
fig7=figure();xlabel('Time (min)'); ylabel('[SciP]');
fig8=figure();xlabel('Time (min)'); ylabel('[CckA_k]');


i=1;
while i<=length(Samples)
    tic
    file=Samples(i);
    baseFileName = matFiles(file).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    matData = load(fullFileName);
    y0=matData.y0;
    try
        Pset=matData.Pbest;
    catch
        Pset=matData.Pset;
    end
    [tout, teout, ieout, yout, ~] = CauloSim('SW','WT', 1600,Pset,y0.');
    [tout,yout,cellarrest] = getLastCycle(tout,yout,teout,ieout);
    if cellarrest == true
        delete(fullFileName);
        i=i+1;
        continue
    end
    K_DivL_CckA = Pset(100); CckAT= yout(:,10); DivL=yout(:,34);DivLDivKP=yout(:,40); DivLDivK=yout(:,19); CckAcdG = yout(:,18);
    DivLT = DivL+DivLDivKP+DivLDivK; CckADivLT = DimerSolve(CckAT,DivLT,K_DivL_CckA); CckADivLDivKP = CckADivLT.*DivLDivKP./DivLT;
    CckAP = CckADivLDivKP+CckAcdG-CckADivLDivKP.*CckAcdG./CckAT; CckAk = CckAT-CckAP;
    
    figure(fig); hold on; plot(tout,yout(:,1)+yout(:,2),'r');
    [~,scalarCtrAT] = findSquares(tout,yout(:,1)+yout(:,2), dpCtrAT);
    CtrATnorm(i)=scalarCtrAT;
    figure(fig1); hold on; 
    plot(tout,yout(:, 31)+yout(:, 18)+2*(yout(:, 38)+yout(:, 39)+yout(:, 41)+yout(:, 43)),'Color',gray)
    figure(fig6); hold on; plot(tout,yout(:, 2),'r');
    [~,scalarCtrAP] = findSquares(tout,yout(:, 2), dpCtrAP);
    CtrAPnorm(i)=scalarCtrAP;
    figure(fig2); hold on; [~,sy] = findSquares(tout,yout(:, 25), dpCpdR); 
    plot(tout,(1/sy)*yout(:, 25),'Color',gray);
    figure(fig3); hold on; [~,sy] = findSquares(tout,yout(:, 3), dpDnaA);
    plot(tout,(1/sy)*yout(:, 3),'b');
    figure(fig4); hold on; [~,sy] = findSquares(tout,yout(:, 4), dpGcrA); plot(tout,(1/sy)*yout(:, 4),'g');
    figure(fig5); hold on; [~,sy1] = findSquares(tout,yout(:, 8), dpCcrM);
    [~,sy2] = findSquares(tout,yout(:, 8), dpCcrM2); plot(tout,((1/sy1+1/sy2)/2)*yout(:, 8),'m');
    figure(fig7); hold on; [~,sy] = findSquares(tout,yout(:, 7), dpSciP); 
    plot(tout,(1/sy)*yout(:, 7),'Color',gold);
    figure(fig8); hold on; [~,sy] = findSquares(tout,CckAk, dpCckAk); 
    plot(tout,(1/sy)*CckAk,'Color',gray);    
    toc
    i=i+1;
end


CtrAPnorm(CtrAPnorm==0)=[]; CtrAPnorm=mean(CtrAPnorm);
CtrATnorm(CtrATnorm==0)=[]; CtrATnorm=mean(CtrATnorm);
figure(fig); hold on; scatter(dpCtrAT(:,1),CtrATnorm*dpCtrAT(:,2), 'ks','filled'); xlim([0 142]);box on;
figure(fig1); hold on; scatter(dpcdG(:,1),dpcdG(:,2), 'ks','filled'); xlim([0 142]);box on;
figure(fig2); hold on; scatter(dpCpdR(:,1),dpCpdR(:,2), 'ks','filled'); xlim([0 142]);box on;
figure(fig3); hold on; scatter(dpDnaA(:,1),dpDnaA(:,2), 'ks','filled'); xlim([0 142]);box on;
figure(fig4); hold on; scatter(dpGcrA(:,1),dpGcrA(:,2), 'ks','filled'); xlim([0 142]);box on;
figure(fig5); hold on; scatter(dpCcrM(:,1),dpCcrM(:,2), 'ks','filled'); xlim([0 142]);box on;
figure(fig5); hold on; scatter(dpCcrM2(:,1),dpCcrM2(:,2), 'k+','linewidth',1.5); xlim([0 142]);box on;
figure(fig6); hold on; scatter(dpCtrAP(:,1),CtrAPnorm*dpCtrAP(:,2), 'ks','filled'); xlim([0 142]);box on;
figure(fig7); hold on; scatter(dpSciP(:,1),dpSciP(:,2), 'ks','filled'); xlim([0 142]);box on;
figure(fig8); hold on; scatter(dpCckAk(:,1),dpCckAk(:,2), 'ks','filled'); xlim([0 142]);box on;

figure(fig)
hold on
xlim([0 142]);
print(fig,['Figures/',File,'A.png'],'-dpng','-r300');  


figure
set(gcf, 'Position', [260 170 453   127*4]);
box on
h(1)=subplot(4,2,1); xlabel('Time (min)');ylabel('[CtrA~P] (\muM)'); xlim([0 142]);
fs=get(gca, 'fontsize'); set(gca, 'fontsize',fs*1.2); box on
h(2)=subplot(4,2,6); xlabel('Time (min)'); ylabel('[cdG] (\muM)'); xlim([0 142]);
fs=get(gca, 'fontsize'); set(gca, 'fontsize',fs*1.2); box on
h(3)=subplot(4,2,7); xlabel('Time (min)'); ylabel('[CpdR]'); xlim([0 142]);
fs=get(gca, 'fontsize'); set(gca, 'fontsize',fs*1.2); box on
h(4)=subplot(4,2,2); xlabel('Time (min)'); ylabel('[DnaA]'); xlim([0 142]); ylim([0 1.2]) 
fs=get(gca, 'fontsize'); set(gca, 'fontsize',fs*1.2); box on
h(5)=subplot(4,2,3); xlabel('Time (min)'); ylabel('[GcrA]'); xlim([0 142]); ylim([0 1.3])
fs=get(gca, 'fontsize'); set(gca, 'fontsize',fs*1.2); box on
h(6)=subplot(4,2,4); xlabel('Time (min)'); ylabel('[CcrM]'); xlim([0 142]); ylim([0 2]);
fs=get(gca, 'fontsize'); set(gca, 'fontsize',fs*1.2); box on
h(7)=subplot(4,2,5); xlabel('Time (min)'); ylabel('[SciP]'); xlim([0 142]);
fs=get(gca, 'fontsize'); set(gca, 'fontsize',fs*1.2); box on
h(8)=subplot(4,2,8); xlabel('Time (min)'); ylabel('[CckA_k]'); xlim([0 142]);
fs=get(gca, 'fontsize'); set(gca, 'fontsize',fs*1.2); box on

% Paste figures on the subplots
copyobj(allchild(get(fig6,'CurrentAxes')),h(1)); box on
copyobj(allchild(get(fig1,'CurrentAxes')),h(2)); box on
copyobj(allchild(get(fig2,'CurrentAxes')),h(3)); box on
copyobj(allchild(get(fig3,'CurrentAxes')),h(4)); box on;
copyobj(allchild(get(fig4,'CurrentAxes')),h(5)); box on
copyobj(allchild(get(fig5,'CurrentAxes')),h(6)); box on
copyobj(allchild(get(fig7,'CurrentAxes')),h(7)); box on
copyobj(allchild(get(fig8,'CurrentAxes')),h(8)); box on
print(gcf,['Figures/',File,'B.png'],'-dpng','-r300');  


end

