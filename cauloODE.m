function dydt = cauloODE(t,y,param,Mutant)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dydt = zeros(49, 1);  

CtrA=y(1); 
CtrAP=y(2); 
DnaAtot=y(3); 
GcrA=y(4);
DivK=y(5);
DivKP=y(6);
SciP=y(7); %y(7) is nothing
CcrM=y(8);
DivJ=y(9);
CckAT=y(10);
PleC=y(11);
PleCk=y(12);
Elong=y(13);
h_cori=y(14);
h_ctra=y(15);
h_ccrm=y(16);
h_perp=y(24);
h_fts=y(45);

V=y(17);
CckAcdG=y(18);
DivLDivK=y(19);
PleCpole=y(21);
PleCtot= y(11)+y(12);
PleCfree=PleCtot-PleCpole;
PodJ=y(22);
PerP=y(23);
CpdR=y(25);
RcdA=y(26);
PopA=y(27);
PdeA=y(28);
PleD=y(29);
PleDP=y(30);
cdG=y(31);
RepSwitch=y(32);
DNAini= y(33);
DivL=y(34);
%SpmX= y(37);
DivJDivK= y(35);
DivJDivKP= y(36);
CpdRP=y(37);
PleDcdG2=y(38);
PleDPcdG2=y(39);
DivLDivKP=y(40);
PopAcdG2= y(41);
X=y(20);
MipZswitch=y(42);
DgcBcdG=y(43);
FtsA= y(44);
FtsAdeg= y(46);
Zring = y(47);
DnaAatp = y(48);
ChromosomeCount = y(49); dydt(49)=0;
% h_FtsZ = y(45);
% DgcB=(39);
% DgcBPde=(40);
% DgcBcdG
% dDgcBdt = -kxcdg*DgCB*

Lon=1;
%kelong = 4*0.95/160*1.22; Pelong=0.05;
kelong = 4*0.95/160*1.22; Pelong=0.05;


ksCtrAP1 = param(1); ksCtrAP2 =  param(2);
kdCtrA1 =  param(3);  kdCtrA2 =  param(4);
JiCtrACtrA =  param(5); JaCtrAGcrA= param(6); JaCtrACtrA=param(7);
JiCtrASciP= param(8); kphosCtrA= param(9);  kdephosCtrA= param(10);
JclpXP=param(11);
epsilonCtrASciP=param(28);
%DnaA
ksDnaA1 = param(12); kdDnaA = param(13);
JiDnaAGcrA = param(14); JiDnaADnaA = param(15);
 JdDnaA=param(17); 
kdDnaAatp= param(114);
%GcrA
ksGcrA = param(18);	kdGcrA = param(19); 
JiGcrACtrA = param(20); JaGcrADnaA=param(21);

%DivK
ksDivK2 =param(22) ; ksDivK= param(23);
kdDivK = param(24); 
JaDivKCtrA = param(25); 
JDivKPPleC = param(26); JDivKPleCk= param(27); 
kdephosDivKP= param(29); kphosDivK=param(30); kphosDivK2=param(31); %kphosDivK2=0.1/2; 


%I and CcrM
JaICtrA = param(34); %%%%***** changed, not changed in spread sheet
ksCcrM = param(35);	kdCcrM = param(36); %%%%***** changed, not changed in spread sheet

%PodJ
ksPodJ=param(37); JPodJGcrA=param(38); JiPodJDnaA= param(39); kd1PodJ=param(40); kd2PodJ=param(41); JdPodJ=param(42);

%PerP
ksPerP=param(43); JaPerPCtrA=param(44); kdPerP=param(45);
% ksPerP=0;F

%PopA
ksPopA=param(46); JiPopAGcrA=param(47); kdPopA=param(48); 
kb_PopAcdg=param(49);

%RcdA
ksRcdA=param(50); kdRcdA=param(51); JaRcdACtrA=param(52); kdRcdA2=param(53); %kdRcdA=0.001;
Krcdapopa=param(54); KclpxpRcdA=param(55);


%PleD
kphosPleDDivJ=param(56); kdephos=param(60);
ksPleD=param(61); JaPleDCtrA=param(62); kdPleD=param(63);

%SciP
ksSciP=param(57); kdSciP=param(58); JaSciPCtrA=param(59);

%PdeA
ksPdeA=param(64); JaPdeACtrA=param(65); kdPdeA1=param(66); kdPdeA2=param(67); JdPdeA=param(68);

%cdG
kscdG=param(69); kdcdG=param(70); Jdcdg=param(71);
kb_xcdg=param(72); ku_xcdg=param(73);
Pde2=param(74);

%PleC
ksPleC = param(75);
ktransPleCk=param(76); ktransPleC=param(77);
kdPleC=param(78); %kdPleCfree=param(79);
%kdPleC=(kdPleCpole*PleCpole/PleCtot+kdPleCfree*PleCfree/PleCtot);
JplecPodJ= param(80);
kPleCbinding= param(81);%
kPleCunbinding=param(82);


%DivJ
ksDivJ = param(83); kdDivJ = param(84);	
ksDivJDivK=param(85); kdDivJDivK=param(86); 
ksDivJDivKP=param(87); kdDivJDivKP=param(88); 

%CpdR
ksCpdR=param(89); kdCpdR=param(90); JaCpdRCtrA=param(91);
kphosCpdR=param(92);  kdephosCpdR=param(93);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DivL
ksDivL=param(94); kdDivL=param(95);
ksDivLDivKP=param(96); kdDivLDivKP=param(97);
ksDivLDivK=param(103);
%CckA
ksCckA=param(98); kdCckA=param(99);
K_DivL_CckA = param(100);
kfcckacdg=param(101); krcckacdg=param(102); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FtsA and Z-ring
k_Zconstrict=param(106);
Jzring=param(107);
ksFtsA=param(108); JaFtsACtrA=param(109); kdFtsA=param(110); 
ksmRNAftsa = param(32); kdmRNA = param(33);

% ksx=param(111); n=cparam(112); c=param(113); Xcap=param(114);
%Cori
 J_coriDnaA=param(116);

KclpxpCpdR=param(117);
DgcB=param(118); kscdG2=param(119);
JiIDnaA = param(120);
kb_xcdg2 = param(122);
% CckAtot=1.3;
% k=35;
JaCpdRDnaA = param(123); JiCpdRGcrA= param(124);
% 
SpmXmutant = param(121);

meth = 1;
methDnaA=param(125); thetaGcrA= param(126); fracCpdR=param(127); methCori = param(128); methftsA=param(129);
thetaCtrA=param(130); methCtrA=param(136); impairfrac=param(137); RiboActrAD51E =param(138);
CtrAbindingmultiplier=param(131); thetaDivJ=param(132); ksCcrMxyl=param(133); mu=param(134); RiboActrA=param(135); 

epsilonPodJgcrA = param(139); epsilonPopA=param(140); methPleC=param(141);
thetaCtrACtrA=param(142); CtrAD51E=param(115);

% %DivK
% ksDivK2 = 0.0028, kdDivK = 0.014, JaDivKCtrA = 2.25
% JDivKPPleC = 0.075, JDivKPleCk=0.15, JDivKDivJ = 0.075
% kdephosDivKP= 17.2, kphosDivK=2.75, kphosDivK2=0.77
% 
% ksDivJ = 0.0025, kdDivJ = 0.025, ksDivJDivK= 200, kdDivJDivK=10
% ksDivJDivKP= 80, kdDivJDivKP=10
% ksPleC = 0.0023, ktransPleCk=0.2, ktransPleC=3200, kdPleC=0.0373
% ksDivL=0.0041; kdDivL=0.0827;
% ksDivLDivKP=600; kdDivLDivKP=10;
% ksCckA=0.009; kdCckA=0.1;
% 
% K_DivL_CckA=1000;
mysteryKinase=param(16);
kscdGsynthetic=0; kdcdgsynthetic=0;

Xcap=(1-thetaDivJ)*(1-SpmXmutant );
DivJactive=(DivJDivK+DivJDivKP)*(X+thetaDivJ)+mysteryKinase;

%% Mutants
FtsAxyl=0; 

%%
% Things Changed:
% CtrA equation, GcrA no longer hill function on P1
% Volume influence changed on PleC localization by PodJ
% Equations for DivJ and PleC phos/kin activity on DivK changed
% got rid of GcrA hill function on PopA 
%got rid of GcrA hill function on CtrA
% added parameter (meth=param(125)) which regulates the degree by which
% methylation impacts gene expression.
%DnaA no longer influences its own promoter
%Addition of DnaA-ATP

%CckA partitions:
DivLT = DivL+DivLDivKP+DivLDivK;
b=-CckAT-DivLT-1/K_DivL_CckA; c=CckAT*DivLT;
CckADivLT= (-b-(b^2-4*c)^(1/2))/2;
% CckADivLT = DimerSolve(CckAT,DivLT,K_DivL_CckA);
CckADivLDivKP = CckADivLT*DivLDivKP/DivLT;
CckAP = CckADivLDivKP+CckAcdG-CckADivLDivKP*CckAcdG/CckAT;
CckAK2 = CckAT-CckAP;
CckAK1 = 0;
% CckAK2 = (DivLDivK/DivLT)*CckADivLT*(1-CckAcdG/CckAT);
% CckAK1 = CckAT-CckAK2-CckAP;

CckAP=CckAP+impairfrac*CckAK1+impairfrac*CckAK2;
CckAK1=(1-impairfrac)*CckAK1;
CckAK2=(1-impairfrac)*CckAK2;

b=-(PopA+PopAcdG2)-RcdA-1/Krcdapopa; c=(PopA+PopAcdG2)*RcdA;
RcdAPopA= (-b-(b^2-4*c)^(1/2))/2;
% RcdAPopA=DimerSolve(PopA+PopAcdG2, RcdA, Krcdapopa);


ClpXP= (CpdR/(CpdR+KclpxpCpdR/V))*RcdAPopA*PopAcdG2/(PopAcdG2+PopA);

% if CtrASciPregRatio > 1 || CtrASciPregRatio < 0
%     error('invalid CtrASciPregRatio')
% end

% ctraP1ratio=0.8;
% dCtrAdt = RiboActrA+(h_ctra*ksCtrAP1*(JiCtrACtrA^2/(JiCtrACtrA^2 + (CtrA+CtrAP)^2))*GcrA/(JaCtrAGcrA + GcrA) + ksCtrAP2*(CtrA+CtrAP)^4/(JaCtrACtrA^4 + (CtrA+CtrAP)^4))*(0.3+0.7*JiCtrASciP^2/(JiCtrASciP^2+SciP^2))- (kdCtrA1+mu)*CtrA - (kdCtrA2*ClpXP*CtrA/(JclpXP+CtrA+CtrAP)) + kdephosCtrA*CckAP*CtrAP/(CtrAP+JCtrACckA) - kphosCtrA*CckAK*CtrA/(CtrA+JCtrACckA);
% dCtrAPdt = -(kdCtrA1+mu)*CtrAP - (kdCtrA2*ClpXP*CtrAP/(JclpXP+CtrA+CtrAP)) - kdephosCtrA*CckAP*CtrAP/(CtrAP+JCtrACckA) + kphosCtrA*CckAK*CtrA/(CtrA+JCtrACckA);
% if CtrA+CtrAP < 3.736e-4
%     CtrAP2=0;
% else
CtrAP2 = ksCtrAP2*(1-thetaCtrACtrA+thetaCtrACtrA*(CtrAP)^2/(JaCtrACtrA^2 + (CtrAP)^2))*JiCtrASciP^2/(JiCtrASciP^2+SciP^2);
% end
CtrAP1 = (1-methCtrA+methCtrA*(2*(1-h_ctra)))*ksCtrAP1*(JiCtrACtrA^2/(JiCtrACtrA^2 + (CtrAP)^2))*(1-thetaCtrA+thetaCtrA*GcrA/(JaCtrAGcrA + GcrA));
dCtrAdt = RiboActrA + CtrAP1+ CtrAP2 - (kdCtrA1+mu)*CtrA - (kdCtrA2*ClpXP*CtrA/(JclpXP+CtrA+CtrAP)) + kdephosCtrA*CckAP*CtrAP - kphosCtrA*CckAK2*CtrA;
dCtrAPdt = -(kdCtrA1+mu)*CtrAP + RiboActrAD51E - (kdCtrA2*ClpXP*CtrAP/(JclpXP+CtrA+CtrAP)) - kdephosCtrA*CckAP*CtrAP + kphosCtrA*CckAK2*CtrA;
% if t> 400
%     dCtrAPdt
%     CtrAP
%     CtrA
%     disp(strcat('part 1 =',num2str(-(kdCtrA1+mu)*CtrAP)))
%     disp(strcat('part 2 =',num2str(- (kdCtrA2*ClpXP*CtrAP/(JclpXP+CtrA+CtrAP)))))
%     disp(strcat('part 3 =',num2str(- kdephosCtrA*CckAP*CtrAP)))
%     disp(strcat('part 4 =',num2str(kphosCtrA*(CckAK1+2.5*CckAK2)*CtrA)))
% end
%JiDnaADnaA/(JiDnaADnaA + DnaA)*

dDnaAtotdt = (ksDnaA1*JiDnaAGcrA/(JiDnaAGcrA + GcrA))*(1-methDnaA+methDnaA*(2*h_cori-1)) ...
     - (kdDnaA+mu)*DnaAtot;
DnaAadp=DnaAtot-DnaAatp;
dDnaAatpdt = (ksDnaA1*JiDnaAGcrA/(JiDnaAGcrA + GcrA))*(1-methDnaA+methDnaA*(2*h_cori-1)) - (kdDnaA+mu+kdDnaAatp*RepSwitch)*DnaAatp;
dGcrAdt = (ksGcrA*JiGcrACtrA^2/(JiGcrACtrA^2 + CtrAP^2))*((1-thetaGcrA)+thetaGcrA*DnaAadp/(DnaAadp+JaGcrADnaA)) - (kdGcrA+mu)*GcrA;
dDivKdt = ksDivK + ksDivK2*CtrAP^2/(JaDivKCtrA^2 + CtrAP^2) + kdephosDivKP*PleC*DivKP ...
    - kphosDivK2*PleCk*DivK - kphosDivK*DivJactive*DivK - (kdDivK+mu)*DivK ...
    - (ksDivJDivK*DivJ*DivK-kdDivJDivK*DivJDivK - kdDivJ*DivJDivK) -(ksDivLDivK*DivL*DivK - kdDivLDivKP*DivLDivK - kdDivL*DivLDivK); 
dDivKPdt = -kdephosDivKP*PleC*DivKP  + kphosDivK2*PleCk*DivK ...
    + kphosDivK*DivJactive*DivK - (kdDivK+mu)*DivKP + 2*(ktransPleCk*PleCk - ktransPleC*PleC*DivKP^2 + kdPleC*PleCk)...
    + 2*kdDivK*PleCk - (ksDivJDivKP*DivJ*DivKP-kdDivJDivKP*DivJDivKP - kdDivJ*DivJDivKP) ...
    -ksDivLDivKP*DivL*DivKP+kdDivLDivKP*DivLDivKP + kdDivL*DivLDivKP; %Since stability of PleC is determined by localization and DivK kd is much smaller, it is better to model release of DivK after PleCk degredation.

dCcrMdt = ksCcrMxyl + ksCcrM*(1-meth+meth*2*(1-h_ccrm))*CtrAP^2/(JaICtrA^2 + CtrAP^2)*JiIDnaA/(JiIDnaA+DnaAtot) - (kdCcrM+mu)*CcrM;
dDivJdt = ksDivJ - (kdDivJ+mu)*DivJ - ksDivJDivK*DivJ*DivK -ksDivJDivKP*DivJ*DivKP + (kdDivJDivKP + kdDivK)*DivJDivKP + (kdDivJDivK + kdDivK)*DivJDivK;
dDivJDivKdt= ksDivJDivK*DivJ*DivK-kdDivJDivK*DivJDivK - (kdDivJ + kdDivK + mu)*DivJDivK;
dDivJDivKPdt= ksDivJDivKP*DivJ*DivKP-kdDivJDivKP*DivJDivKP - (kdDivJ + kdDivK + mu)*DivJDivKP;
dDivLdt= ksDivL - (kdDivL+mu)*DivL -(ksDivLDivK*DivL*DivK - kdDivLDivKP*DivLDivK) - ksDivLDivKP*DivL*DivKP + kdDivLDivKP*DivLDivKP + kdDivK*(DivLDivKP+DivLDivK);
dDivLDivKPdt= ksDivLDivKP*DivL*DivKP-kdDivLDivKP*DivLDivKP - (kdDivL + kdDivK + mu)*DivLDivKP;
% ksDivLDivK=0;fts
dDivLDivKdt = ksDivLDivK*DivL*DivK - kdDivLDivKP*DivLDivK - (kdDivL + kdDivK + mu)*DivLDivK;

dCckATdt = ksCckA - (kdCckA+mu)*CckAT; %- ktransCckAk*CckAk + (ktransCckA*(CckAtot - CckAk)/(CckAtot - CckAk + JCckADivL))*DivLf;


dPleCdt = ksPleC*(1-methPleC+methPleC*(2*(1-h_fts))) - (kdPleC+mu)*PleC + ktransPleCk*PleCk - ktransPleC*PleC*DivKP^2 + 2*kdDivK*PleCk ;
dPleCkdt = -(kdPleC + 2*kdDivK + mu)*PleCk - ktransPleCk*PleCk + ktransPleC*PleC*DivKP^2;

dElongdt = kelong*RepSwitch;
% dh_coridt = - kmccrM*CcrM^10/(JmccrM^10 + CcrM^10)*h_cori;
% dh_ctradt =  - kmccrM*CcrM^10/(JmccrM^10 + CcrM^10)*h_ctra;
% dh_ccrmdt = - kmccrM*CcrM^10/(JmccrM2^10 + CcrM^10)*h_ccrm;
% dh_ftsdt =  - kmccrM*CcrM^10/(JmccrM^10 + CcrM^10)*h_fts;
% dh_PerPdt = -kmccrM*CcrM^10/(JmccrM^10 + CcrM^10)*h_perp;

% M=1-hcori*0.5
% hcori=2*(1-M)
dh_coridt = 0;
dh_ctradt = 0;
dh_ccrmdt = 0;
dh_ftsdt =  0;
dh_PerPdt = 0;


dVdt = mu*V;
dZdt=0;

% dPleCpoledt= kPleCbinding/V*PleCfree*PodJ^2/(PodJ^2+(JplecPodJ/V)^2)-kPleCunbinding*PleCSig*PleCpole-(kdPleCpole+mu)*PleCpole;
dPleCpoledt= kPleCbinding*PleCfree*PodJ/(PodJ*V+JplecPodJ)-kPleCunbinding*PleCpole-(kdPleC+mu)*PleCpole;


dPodJdt= ksPodJ*(1-meth+meth*(2*(1-h_perp)))*(1-epsilonPodJgcrA+epsilonPodJgcrA*GcrA/(JPodJGcrA+GcrA))*(JiPodJDnaA/(DnaAtot+JiPodJDnaA))-(kd1PodJ+mu)*PodJ-kd2PodJ*PodJ*PerP/(PodJ+JdPodJ);
dPerPdt= (ksPerP*CtrAP^2/(CtrAP^2+JaPerPCtrA))*(1-meth+meth*(2*(1-h_perp))) - (kdPerP+mu)*PerP;
dCpdRdt= ksCpdR*CtrAP/(JaCpdRCtrA+CtrAP)*((1-fracCpdR)+fracCpdR*DnaAtot/(JaCpdRDnaA+DnaAtot))*(JiCpdRGcrA/(GcrA + JiCpdRGcrA)) + kdephosCpdR*CckAP*CpdRP - kphosCpdR*CckAK2*CpdR - (kdCpdR+mu)*CpdR;
dCpdRPdt= -kdephosCpdR*CckAP*CpdRP + kphosCpdR*CckAK2*CpdR - (kdCpdR+mu)*CpdRP;
dRcdAdt= ksRcdA*CtrAP^2/(CtrAP^2+JaRcdACtrA^2)-(kdRcdA+mu)*RcdA - kdRcdA2*CpdR*(RcdA-RcdAPopA)/(CpdR+KclpxpRcdA);

dPopAdt= ksPopA*(1-epsilonPopA+epsilonPopA*JiPopAGcrA/(JiPopAGcrA+GcrA))-(kdPopA+mu)*PopA - kb_PopAcdg*PopA*cdG^2 + ku_xcdg*PopAcdG2;
dPopAcdG2dt= kb_PopAcdg*PopA*cdG^2 - ku_xcdg*PopAcdG2 - (kdPopA+mu)*PopAcdG2;

dPdeAdt= ksPdeA*CtrAP/(CtrAP+JaPdeACtrA) - (kdPdeA1+mu)*PdeA - kdPdeA2*CpdR*PdeA/(JdPdeA + PdeA); 

dPleDdt= ksPleD*CtrAP^2/(CtrAP^2+JaPleDCtrA^2) - (kdPleD+mu)*PleD - kphosPleDDivJ*DivJactive*PleD ...
    + ((1/10)*kdephos*PleC+ (9/10)*kdephos*PleC*PleCfree/PleCtot)*PleDP ...
    - kb_xcdg*PleD*cdG^2 + ku_xcdg*PleDcdG2;

dPleDPdt= kphosPleDDivJ*DivJactive*PleD - ((1/10)*kdephos*PleC+ (9/10)*kdephos*PleC*PleCfree/PleCtot)*PleDP  ...
    - (kdPleD+mu)*PleDP - kb_xcdg*PleDP*cdG^2 + ku_xcdg*PleDPcdG2; 
dPleDcdG2dt= kb_xcdg*PleD*cdG^2 - ku_xcdg*PleDcdG2 - (kdPleD+mu)*PleDcdG2 - kphosPleDDivJ*DivJactive*PleDcdG2 ...
     + ((1/10)*kdephos*PleC+ (9/10)*kdephos*PleC*PleCfree/PleCtot)*PleDPcdG2 ;
dPleDPcdG2dt= kb_xcdg*PleDP*cdG^2 - ku_xcdg*PleDPcdG2 - (kdPleD+mu)*PleDPcdG2 + kphosPleDDivJ*DivJactive*PleDcdG2 ...
     - ((1/10)*kdephos*PleC+ (9/10)*kdephos*PleC*PleCfree/PleCtot)*PleDPcdG2 ;

%Changed equation... added DgcB as another DGC enzyme. DgcB is inhibited by
%PdeA... assumed tight binding. 
dcdGdt= kscdGsynthetic-kdcdgsynthetic*cdG+kscdG*PleDP+kscdG2*max(DgcB-PdeA,0)*((DgcB-DgcBcdG)/DgcB)- kdcdG*(PdeA+Pde2)*cdG/(cdG+Jdcdg) -mu*cdG - 2*(kb_PopAcdg*PopA*cdG^2 - ku_xcdg*PopAcdG2 - kdPopA*PopAcdG2) + ... 
    2*(-kb_xcdg*(PleD+PleDP)*cdG^2+(kdPleD+ku_xcdg)*(PleDcdG2+PleDPcdG2)) + (-kfcckacdg*cdG*(CckAT-CckAcdG) + krcckacdg*CckAcdG + kdCckA*CckAcdG) ...
    - 2*kb_xcdg2*(DgcB-DgcBcdG)*cdG^2+2*ku_xcdg*DgcBcdG;

dCckAcdGdt = kfcckacdg*cdG*(CckAT-CckAcdG) - krcckacdg*CckAcdG - (kdCckA+mu)*CckAcdG;
dDgcBcdGdt = kb_xcdg2*(DgcB-DgcBcdG)*cdG^2-ku_xcdg*DgcBcdG;

%X
ksx=0.14; n=3; c=0.035;
% dXdt= ksx*(X^n/(X^n+c^n)-X^n/(Xcap^n+c^n));
dXdt= ksx*(X^n/(X^n+c^n))*(Xcap-X);
% if Elong < 0.025

% dDNAinidt= 0.15*(1-0.5+0.5*(2*h_cori-1))*(J_coriCtrA^5/(J_coriCtrA^5+CtrAP^5))*(DnaAatp^4/(J_coriDnaA^4+DnaAatp^4))-0.02*DNAini;
% try
%     [DNAf, DNACtrA,DNACtrAP,DNA2CtrA2,DNA2CtrAP2,DNA2CtrACtrAP] = DNAbinding(CtrAbindingmultiplier*CtrA*10,CtrAP*10,0.979,0.0337,2.834e-5);
% catch
%     Mutant
%     t
%     [DNAf, DNACtrA,DNACtrAP,DNA2CtrA2,DNA2CtrAP2,DNA2CtrACtrAP] = DNAbinding(CtrAbindingmultiplier*CtrA*10,CtrAP*10,0.979,0.0337,2.834e-5);
% end

% k=0.00806;
% k2=0.0366;
% DNA2CtrAP2=  CtrAP^2/(k^2+CtrAP^2+k2^2*CtrA^2);
kd1=1.065; kd2=0.039; kd3=0.000085*(1-CtrAD51E)+kd2*CtrAD51E;
DNAF=kd1*(kd1+CtrAbindingmultiplier*2*CtrA+2*CtrAP+(CtrAbindingmultiplier*CtrA)^2/kd2+CtrAP^2/kd3+CtrAbindingmultiplier*2*CtrA*CtrAP/kd2)^-1;
DNA2CtrAP2= (DNAF*CtrAP^2)/(kd1*kd3);
dDNAinidt= (1-methCori+methCori*(2*h_cori-1))*(1-DNA2CtrAP2)^5*(DnaAatp/(J_coriDnaA+DnaAatp))^2-0.05*DNAini;

% else 
%     dDNAinidt=0;
% end
% dDNAinidt= 10*(1-0.5+0.5*(2*h_cori-1))*(J_coriCtrA^5/(J_coriCtrA^5+CtrAP^5))*(DnaAatp^4/(J_coriDnaA^4+DnaAatp^4))-0.1*DNAini;
dRepSwitchdt=0;


%%Z-ring protiens

dMipZswitchdt=0;
dPleCSigdt =  0;
dFtsAdt = ksFtsA*(CtrAP^2/(CtrAP^2+JaFtsACtrA^2))*(1-methftsA+methftsA*2*(1-h_fts))-(kdFtsA+mu+FtsAdeg)*FtsA +FtsAxyl; 

dZringdt= -k_Zconstrict*MipZswitch*FtsA^5/(FtsA^5+(Jzring+3*Zring)^5);

dSciPdt= ksSciP*(CtrAP^2/((CtrAP)^2+JaSciPCtrA^2))-(kdSciP+mu)*SciP;
dFtsAdegdt=0;

dydt(1)=dCtrAdt; dydt(2)=dCtrAPdt; dydt(3)=dDnaAtotdt; dydt(4)=dGcrAdt;
dydt(5)=dDivKdt; dydt(6)=dDivKPdt; dydt(7)=dSciPdt; dydt(8)=dCcrMdt;
dydt(9)=dDivJdt; dydt(10)=dCckATdt; dydt(11)=dPleCdt; dydt(12)=dPleCkdt;
dydt(13)=dElongdt; dydt(14)=dh_coridt; dydt(15)=dh_ctradt; 
dydt(16)=dh_ccrmdt; dydt(17)=dVdt; dydt(18)=dCckAcdGdt; dydt(19)=dDivLDivKdt;
dydt(20)=dXdt; dydt(21)=dPleCpoledt; dydt(22)=dPodJdt; dydt(23)=dPerPdt;
dydt(24)=dh_PerPdt; dydt(25)=dCpdRdt; dydt(26)=dRcdAdt; dydt(27)=dPopAdt; 
dydt(28)=dPdeAdt; dydt(29)=dPleDdt; dydt(30)=dPleDPdt; dydt(31)=dcdGdt;
dydt(35)=dDivJDivKdt; dydt(36)=dDivJDivKPdt; dydt(37)=dCpdRPdt; 
dydt(38) = dPleDcdG2dt; dydt(39) = dPleDPcdG2dt; dydt(40)=dDivLDivKPdt;
dydt(41)= dPopAcdG2dt; dydt(42)= dMipZswitchdt; dydt(43)= dDgcBcdGdt; dydt(44)= dFtsAdt;
dydt(45)=dh_ftsdt; dydt(46)=dFtsAdegdt; dydt(47)= dZringdt;  dydt(48)=dDnaAatpdt;
dydt(32)=dRepSwitchdt;
dydt(33)=dDNAinidt;
dydt(34)=dDivLdt;

end

