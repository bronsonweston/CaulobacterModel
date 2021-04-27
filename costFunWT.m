function [cost] = costFunWT(t, y, param,CtrAversion)
%Calculates cost of simulation

scaleCtrAP=2; scaleCtrAT=100; scaleDnaA=10; scaleGcrA=10;
scalePleCpole=5; scaleCcrM=15;scalePodJ=3;scalePerP=3;
scaleCckAK=3;scaleCdG=1000;scaleRcdA=2;scaleCpdR=5; scaleCpdRT=7.5;
scaleSciP=45; scaleFtsA=10; scalePdeA=10; scaleDivJ=1;
scaleCCT=1; scaleCRT=1; scaleDivKP=10; scalePleC=30; scalePleCtot=200;


K_DivL_CckA = param(100); KclpxpCpdR = param(117); Krcdapopa=param(54);
CckAT= y(:,10); DivL=y(:,34);DivLDivKP=y(:,40); DivLDivK=y(:,19); 
CckAcdG = y(:,18); PopA= y(:,27); PopAcdG2=y(:,41);
RcdA = y(:,26); CpdR= y(:,25); V= y(:,17);

DivLT = DivL+DivLDivKP+DivLDivK;
b=-CckAT-DivLT-1./K_DivL_CckA; c=CckAT.*DivLT;
CckADivLT= (-b-(b.^2-4*c).^(1/2))./2;
% CckADivLT = DimerSolve(CckAT,DivLT,K_DivL_CckA);
CckADivLDivKP = CckADivLT.*DivLDivKP./DivLT;
CckAP = CckADivLDivKP+CckAcdG-CckADivLDivKP.*CckAcdG./CckAT;
CckAK2 = (DivLDivK./DivLT).*CckADivLT.*(1-CckAcdG./CckAT);
CckAK1 = CckAT-CckAK2-CckAP;
CckAK=CckAK1+CckAK2;
RcdAPopA=DimerSolve(PopA+PopAcdG2, RcdA, Krcdapopa);
ClpXP= (CpdR./(CpdR+KclpxpCpdR./V)).*RcdAPopA.*PopAcdG2./(PopAcdG2+PopA);

%% Data Points
dpCtrAP = [10*1.07 0.4; 80*1.07 0.5; 105*1.07 .8; 130*1.07 1]; %Jacobs 2003 (relative to max = 1) removed second point (t=40, P=0.15) as CtrAT suggests it is incorrect
dpCtrAP(:,1)=dpCtrAP(:,1)*t(end-1)/130;
if strcmp(CtrAversion, 'Complete')
    dpCtrAT= [0 0.8; 20 0; 20 0; 20 0; 40 0; 60 0; 80 0.5; 100 1; 120 .85; 140 0.7]; %Collier 2006 - DnaA couples DNA... (Normalized based on max expression)
    dpCtrAT(:,1)=dpCtrAT(:,1)*t(end-1)/140;
    dpCtrAT2= [ 18 0; 20 0; 20 0;20 1e-3];
    damp=1/100;
elseif strcmp(CtrAversion, 'Slow') || strcmp(CtrAversion, 'CtrAbindingOff')
    Mcgrath = [0	12915.912; 20	9480.648; 40	2687.648; 60	515.698; 80	5769.719; 100	10582.426; 120	13602.205;140	13106.439];
    Mcgrath(:,2)=Mcgrath(:,2)/max(Mcgrath(:,2));
    dpCtrAT = Mcgrath;
    dpCtrAT(:,1)=dpCtrAT(:,1)*t(end-1)/140;
    dpCtrAT2= [40 4.5; 60 0; 60 1e-3];
%         dpCtrAT = [0 0.6; 20 0.2; 20 0.2; 40 0; 60 0.05; 80 0.45; 100 .95; 120 1]; % Holtzendorff-2004- Oscillating global regulatiors
%     dpCtrAT(:,1)=dpCtrAT(:,1)*t(end-1)/130;
%     dpCtrAT2= [40 0; 40 0; 40 0; 40 0;40 1e-3];
    damp=1/2;
else
    error('Invalid Type of Parameter Collection Specified')
end

dpSciP = [0 0.87; 20 0.856; 40 0.043; 80 0; 100 0.02; 120 0.4; 140 1]; %Tan et al. 2010
dpSciP(:,1)=dpSciP(:,1)*t(end-1)/140;

% dpCtrAT3 = [0 1; 5 0.85; 10 0.55; 20 0.38; 25 0.25; 30 0.1; 40 0.05; 50 0.05; 60 0; 70 0; 80 0.05; 90 0.1; 100 0.2; % Quon et al. 1997
% dpDnaA = [0 .5; 20 1; 40 0.9; 60 0.3; 80 0.15; 100 .1; 120 0.4; 140 0.6]; %Collier 2006 - DnaA couples DNA... (Normalized based on max expression)
dpDnaA = [0 0.78;15 1; 30 0.9; 45 0.725; 60 0.375; 75 0.325; 90 0.225; 105 0.42; 120 0.4; 135 0.6]; % Cheng and Keiler 2009 (cell cycle time 136 +- 11 min)
dpDnaA(:,1) = dpDnaA(:,1)*t(end-1)/136;
dpGcrA2 = [0 0; 20 .5; 40 1; 60  0.9; 80 0.98; 100 0.58; 120 .65; 140 .75];
dpGcrA = [0 0; 20 .2; 40 .7;60  1; 80 0.8; 100 0.6; 120 .5; 140 .4]; % Holtzendorff-2004- Oscillating global regulatiors
dpPleCpole = [12*0.93 0.68; 28*0.93 0.37; 40*0.93 0.25; 60*0.93 0.4; 75*0.93 0.5; 90*0.93 0.72; 105*0.93 0.85; 120*0.93 1; 130*0.93 0.92; 150*0.93 0.72]; %Christen 2010
Zcloseset=find(y(:,47)<=0);
tclose=t(Zcloseset(2));
plecpoints=find(dpPleCpole(:,1)>tclose);
dpPleCpole(plecpoints,2)=dpPleCpole(plecpoints,2)./0.46;
dpPleCtot = [0 .85/0.85; 20 0.5/0.85; 40 .28/0.85; 60 .38/0.85; 80 .45/0.85; 100 .65/0.85; 120 0.75/0.85; 140 0.8/0.85]; %Viollier 2002 - ommitted last datapoint bc it was at time point 160 and is unlikely behavior
dpPodJ1 = [0 0; 20 0; 40 .1; 60 .55; 80 1; 100 1; 120 0.42; 140 0.18]; %Chen et al.
dpPodJ2 = [0 0; 20 0.05; 40 0.1; 60 0.9; 80 1; 100 0.85; 120 0.27; 140 0.15]; %Vollier et al.
dpPerP = [0*0.93 0.6; 20*0.93 0.25; 40*0.93 0; 60*0.93 0; 60*0.93 0; 95*0.93 0.5; 120*0.93 1; 130*0.93 0.8; 150*0.93 0.6]; % Shengua Li's model
dpPerP(:,1)= dpPerP(:,1)*t(end-1)/150;
dpDivJ = [0 0.364; 20 0.574; 40 0.81; 60 0.997; 80 0.924; 100 0.9; 120 1; 140 0.94]; % Wheeler et al. 1999
dpDivJ2 = [0 0.26; 20 0.56; 40 0.86; 60 0.81; 80 0.97; 100 0.93; 120 0.86; 140 1]; %Sanselicio 2015
dpDivKP = [0 0.68; 30 0.87; 60 0.87; 90 0.96; 120 1]; %Jacobs et al. 2001 (skewed from a 100 min cell cycle)
dpDivKtot = [0 0.88; 25 0.84; 50 0.79; 75 0.98; 100 0.98; 125 1]; %Jacobs et al. 2001 (skewed from a 120 min cell cycle)
%Rough estimate from Hubert Lam 2003 (0 means low or min level)
dpCckAk = [10 0.4; 40 0.3; 70 0.4; 100 1; 130 1]; %Jacobs 2003 (relative to max = 1)
dpCcrM = [5 .15; 45 0; 90 0.5; 120 1; 120 1; 150 0.25]; % Grunenfelder et al. 2001 Proteomic analysis of the bacterial cell cycle (relative to max)
dpPerP(:,1)= dpPerP(:,1)*t(end-1)/150;
dpDNA = [0 0; 5 0; 10 0; 15 0.05; 20 0.18; 25 0.25; 30 0.4; 40 0.65; 50 0.75; 60 0.9; 70 0.85; 80 1; 90 0.9; 100 0.95]; % Quon et al. 1997
dpCpdRtot = [1*0.93 0.97; 25*0.93 1; 50*0.93 0.66; 75*0.93 0.497; 100*0.93 0.502; 125*0.93 0.64; 150*0.93 0.96]; %Iniesta et al. 2006 Figure 5A
dpCpdRtot(:,1)= dpCpdRtot(:,1)*t(end-1)/150;
dpCpdR = [1*0.93 0.42; 25*0.93 0.938; 50*0.93 0.636; 75*0.93 0.48; 100*0.93 0.0; 105 0; 150*0.93 0]; %Iniesta et al. 2006 Figure 5A estimated
dpCpdR(:,1)= dpCpdR(:,1)*t(end-1)/140;
dpcdG = [0 1/3; 20 1; 20 1; 20 1; 40 0.55; 40 0.55; 40 0.55; 60 0.35; 80 0.3; 100 0.25; 120 0.05; 140 0.35]; %Abel 2013
dpcdG(:,1)= dpcdG(:,1)*t(end-1)/140;
dpRcdA=[0 0.1; 20 0.95; 40 1; 60 0.55; 80 0.55; 100 0.5; 120 0.4; 140 0.7]; %Mcgrath 2006
dpRcdA(:,1)= dpRcdA(:,1)*t(end-1)/140;
dpPdeA = [0 1; 18.75 0.41; 37.5 0.158; 56.25 0.065; 75 0.072; 93.75 0.173; 112.5 0.44; 131.25 0.82; 150 0.998]; %Abel 2011 (my analysis of western)
dpPdeA(:,1)= dpPdeA(:,1)*t(end-1)/150;
dpFtsA = [0 250; 15 150; 30 0; 45 0; 60 355; 75 450; 90 800; 105 1700; 120 2300; 135 3200; 150 1600]; %in molecules Martin 2004
dpFtsA(:,1)= dpFtsA(:,1)*t(end-1)/150;

% dpSciP = [0 0.87; 20 0.856; 40 0.043; 80 0; 100 0.02; 120 0.4; 140 1]; %Tan et al. 2010
imaginedpPleDP = [25 1; 100 0];
imaginedpDivKP = [25 1; 100 0.1];
imaginedpPleCPhosFrac = [25 0; 120 0.8824];

DivLTsum=trapz(t,DivLT)/t(end);
DivLfreesum=trapz(t,DivL)/t(end);
divLratio=DivLfreesum/DivLTsum;
if divLratio > 0.5
    DivLcost= ((divLratio-0.5)*10)^2;
else
    DivLcost=0;
end

%dpPleD =
%dpPleCP = [15 1; 45 0; 80 0; 135 1];


FtsAc=scaleFtsA*findSquares(t,y(:,44),dpFtsA);
PdeAc=scalePdeA*findSquares(t,y(:,28),dpPdeA,0.3);
RcdAc=scaleRcdA*findSquares(t,y(:,26),dpRcdA,0.15);
CckAKc=scaleCckAK*findSquares(t,CckAK,dpCckAk);
CpdRTc = scaleCpdRT*findSquares(t,y(:, 25)+y(:, 37),dpCpdRtot,1);
CpdRc =scaleCpdR*findSquares(t,y(:, 25),dpCpdR);
CdGc=scaleCdG*findSquares(t,y(:, 31)+y(:, 18)+2*(y(:, 38)+y(:, 39)+y(:, 41)),dpcdG,0.27);
CtrAPc=scaleCtrAP*findSquares(t,y(:,2),dpCtrAP);
CtrATc=scaleCtrAT*findSquares(t,y(:,1)+y(:,2),dpCtrAT,25);
CtrAT2c=damp*findSquares2(t,y(:,1)+y(:,2),dpCtrAT2);
DnaAc =scaleDnaA*findSquares(t,y(:,3),dpDnaA,2);
GcrAc =scaleGcrA*(findSquares(t,y(:,4),dpGcrA,1.5)+findSquares(t,y(:,4),dpGcrA2,1.5))/2;
[PleCtotc,PleCtotScale] = findSquares(t,y(:,11)+y(:,12),dpPleCtot);
PleCtotc=scalePleCtot*PleCtotc;
PleCpolec = scalePleCpole*findSquares(t,y(:,21),dpPleCpole,0.8*PleCtotScale);
CcrMc =scaleCcrM*findSquares(t,y(:,8),dpCcrM);
SciPc =scaleSciP*findSquares(t,y(:,7),dpSciP,1.25);
% SciPc=0;
PodJc =scalePodJ*(findSquares(t,y(:,22),dpPodJ1,1)+findSquares(t,y(:,22),dpPodJ2,1))/2;
PerPc =scalePerP*findSquares(t,y(:,23),dpPerP,1);
DivJc =scaleDivJ*(findSquares(t,y(:,9),dpDivJ,0.05)+findSquares(t,y(:,9),dpDivJ2,0.05))/2;
PleCPhosc = scalePleC*findSquares(t,y(:,11)./(y(:,11)+y(:,12)),imaginedpPleCPhosFrac,1);
% DivKtotscale = max(y(:, 5)+y(:, 6)+2*y(:, 12)+y(:, 35)+y(:, 36)+y(:, 40));
DivKPc = scaleDivKP*findSquares(t,(y(:, 6)+2*y(:, 12)+y(:, 36)+y(:, 40))./(y(:, 5)+y(:,19)+y(:, 6)+2*y(:, 12)+y(:, 35)+y(:, 36)+y(:, 40)),imaginedpDivKP,1);
% PleDtotscale = max(y(:, 29)+y(:, 30)+y(:, 38)+y(:, 39));
PleDPc = findSquares(t,y(:, 30)+y(:, 39),imaginedpPleDP,0.025);


% ctraratiopoint=find(y(:,2)==max(y(:,2)));
% CtrA2Pratiocost=(15*(1-CtrA2Pratio))^2
A=find(y(:,47)<=0);
B=find(t>(t(A(2))+5));
t(B(1));
CtrA2Pratio= y(B(1),2)/(y(B(1),1)+y(B(1),2));
CtrA2Pratiocost=(10*(1-CtrA2Pratio))^2;
% FtsAc
% PdeAc
% RcdAc
% CckAKc
% CpdRTc
% CpdRc
% CdGc
% CtrAPc
% CtrATc
% DnaAc
% GcrAc
% PleCpolec
% PleCtotc 
% CcrMc
% PodJc
% PerPc
% DivJc
% PleCPhosc
% DivKPc
% PleDPc
% DivLcost
% CtrAT2c
concentrationCosts= FtsAc+ PdeAc+ RcdAc+ CckAKc+ CpdRTc+ CtrAT2c+...
    CpdRc+ CdGc+ CtrAPc+ CtrATc+ DnaAc + GcrAc + PleCpolec + PleCtotc + SciPc+ ...
    CcrMc + PodJc + PerPc + DivJc +PleCPhosc+DivKPc+PleDPc + CtrA2Pratiocost;% + DivLcost;

ChromosomeRepTime=min(t(y(:,13)>0));
CellCycleTime = t(end);
% ChromosomeRepTime
% scaleCRT*((ChromosomeRepTime-20)/2)^2 
% CellCycleTime
% scaleCCT*((CellCycleTime-140)/8)^2
pks=findpeaks(y(:, 33));
points=ismember(y(:,33),pks);
pksindex = find(points);
if length(pksindex)>1
    if pksindex(2) >0.35
        IniCost2= ((pks(2)-0.35)*10)^2;
    else
        Inicost2=0;
    end
end
if length(pksindex)>2    
    IniCost3= (pks(3)*10)^2;
else
    IniCost3=0;
end
eventCosts= scaleCRT*((ChromosomeRepTime-18)/0.8)^2 + scaleCCT*((CellCycleTime-145)/3)^2;

cost= concentrationCosts + eventCosts +  IniCost2 + IniCost3;
end

