function [param] = getMutantParam(Mutant,Pset)
param=Pset;
if strcmp(Mutant,'DNAbindingCtrAoff')
    param(131)=0;
elseif strcmp(Mutant,'DNAbindingCtrAoff & CtrAdelta3omega')
    param(131)=0;
    param(1)=0; param(2)=0;
    param(135)=25*(Pset(3)+Pset(134));
    param(4)=0;
elseif strcmp(Mutant,'ClearCtrA')
    param(3)=10000;
elseif strcmp(Mutant, 'deltaSpmX')
    param(121)=1;
elseif strcmp(Mutant, 'LS1')
    param(121)=0.9;
    param(133)=param(35);
elseif strcmp(Mutant,'LS1+PpleC::Tn')
    param(75)=0;
    param(121)=0.98;
    param(133)=param(35);
elseif strcmp(Mutant, 'DNAdamage')
    param(1)=0.0002;
    param(2)=0;
    param(106)=0;
% elseif strcmp(Mutant, 'StarveX')
%     ksx=0;
%     mu=0;
% elseif strcmp(Mutant,'StarveLon')
%     Lon=100;
% %     mu=0;
% %     kscdG=0;
% elseif strcmp(Mutant,'StarveCtrA')
%     kdCtrA2=0;
%     mu=0;
elseif strcmp(Mutant,'deltaGcrA')
    param(18)=0;
elseif strcmp(Mutant,'deltaGcrA+ctrAP1gcrA')
    param(18)=0;
    param(130)=0;
elseif strcmp(Mutant,'deltaGcrA+deltaRcdA')
    param(18)=0;
    param(50)=0;
elseif strcmp(Mutant,'deltaRcdA')
    param(50)=0; param(51)=Pset(51)*100;
elseif strcmp(Mutant,'StarveAll')
    param(121)=1;
    param(12)=0;
    param(69)=0;
    param(134)=0.0054/3;
    param(1)=Pset(1)/3;
    param(2)=Pset(2)/3;
elseif strcmp(Mutant,'StarveAll-cdG0')
    param(121)=1;
    param(12)=0;
    param(4)=0;
    param(134)=0.0054/3;
    param(1)=Pset(1)/2;
    param(2)=Pset(2)/2;
% elseif strcmp(Mutant,'StarveAll+CckA')
%     param(121)=1;
%     ksCtrAP1=ksCtrAP1/2;
%     ksCtrAP2=ksCtrAP2/2;
%     kdCtrA2=0;
%     %     kscdG=0;
% %     mu=0.0054/3;
%     ksCckA=ksCckA*3.5;
% elseif strcmp(Mutant,'StarveAll50%')
% %     ksx=0;
%     Lon=2.4;
%     kdCtrA2=kdCtrA2/4;
%     mu=0.0054/1.5;
%     param(121)=1;
elseif strcmp(Mutant, 'deltaDnaA')
    param(12)=0; param(13)=0;
elseif strcmp(Mutant, 'RiboActrA') 
    param(1)=0; param(2)=0; param(135)=25*(Pset(3)+Pset(134));
elseif strcmp(Mutant, 'CtrAdelta3omega') 
%     param(1)=0; param(2)=0; 
%     param(135)=25*(Pset(3)+Pset(134));
    param(4)=0; 
elseif strcmp(Mutant, 'CtrAD51E')
    param(1)=0; param(2)=0; 
    param(138)=30*(Pset(3)+Pset(134));
    param(10)=0;
    param(115)=1;
elseif strcmp(Mutant, 'CtrAD51E~P')
    param(1)=0; param(2)=0;
    param(138)=30*(Pset(3)+Pset(134));
    param(10)=0;
elseif strcmp(Mutant, 'CtrAD51Edelta3omega') 
    param(1)=0; param(2)=0; 
    param(138)=25*(Pset(3)+Pset(134));
    param(4)=0; param(10)=0;
elseif strcmp(Mutant, 'SM921') 
    param(1)=param(1)/100000;
elseif strcmp(Mutant, '100%GcrAP1') 
    param(130)=1;
elseif strcmp(Mutant, 'deltaGcrA & 100%GcrAP1') 
    param(130)=1;
    param(18)=0;
elseif strcmp(Mutant, 'CtrA_dephos_off')
    param(10)=0;
elseif strcmp(Mutant, 'deltaccrM') 
    param(35)=0;
% elseif strcmp(Mutant, 'deltaccrM & FtsAxyl')
%     param(35)=0;
%     FtsAxyl=ksFtsA/10;
% %     MipZswitch=1;
elseif strcmp(Mutant, 'deltaccrM & cdG0')
    param(35)=0;
    param(69)=0; param(119)=0;
elseif strcmp(Mutant,'cdG0')
    param(69)=0; param(119)=0;
elseif strcmp(Mutant, 'PdivK::Tn & cdG0') || strcmp(Mutant, 'cdG0 & PdivK::Tn')
    param(69)=0; param(119)=0;
    param(22) =Pset(22)*0.1; param(23)=Pset(23)*0.1; 
elseif strcmp(Mutant,'PdivK::Tn')
    param(22) =Pset(22)*0.4; param(23)=Pset(23)*0.4; 
elseif strcmp(Mutant,'PdivK::Tn & divL(Y550F)') || strcmp(Mutant,'divL(Y550F) & PdivK::Tn')
    param(22) =Pset(22)*0.15; param(23)=Pset(23)*0.15;
    %param(96)=Pset(96)*10; param(103)=Pset(103)*1000;
    param(97)=Pset(97)/10;
elseif strcmp(Mutant,'deltapopA & PdivK::Tn') || strcmp(Mutant,'PdivK::Tn & deltapopA')
    param(22) =Pset(22)*0.15; param(23)=Pset(23)*0.15;
    param(46)= 0;
elseif strcmp(Mutant,'deltapopA')
    param(46)= 0;
elseif strcmp(Mutant,'deltapopA & deltaPleD') || strcmp(Mutant,'deltaPleD & deltapopA')
    param(61)=0;
    param(46)= 0;
elseif strcmp(Mutant,'deltapopA & deltaPleD & PdivK::Tn')
    param(22) =Pset(22)*0.25; param(23)=Pset(23)*0.25;
    param(61)=0;
    param(46)= 0;
elseif strcmp(Mutant,'deltaPleD & PdivK::Tn')
    param(22) =Pset(22)*0.25; param(23)=Pset(23)*0.25;
    param(61)=0;
elseif strcmp(Mutant,'deltadgcB')
    param(118)= Pset(118)/100000;
elseif strcmp(Mutant,'PpleC::Tn & deltadivJ')
    param(83)=0;
    param(75)=0;
    param(121)=0.98;
elseif strcmp(Mutant,'deltadivJ')
    param(83)=0;
elseif strcmp(Mutant,'divKxyl')
    param(23)=Pset(23)*25;
elseif strcmp(Mutant,'divKcs')
    param(22) = 0; param(23)=0; param(24)=Pset(24)*50;
elseif strcmp(Mutant,'PpleC::Tn')
    param(75)=0;
    param(121)=0.98;
elseif strcmp(Mutant,'PpleC::Tn & deltaSpmX')
    param(75)=0;
    param(121)=1;
% elseif strcmp(Mutant,'PpleC::Tn + dFtsZ')%read up more on this
%     ksPleC=0;
%     param(121)=0.9;
%     ksFtsA=0;
% elseif strcmp(Mutant,'PpleC::Tn Alternate')%read up more on this
%     ksPleC=0;
% %     Xcap=0.1;
elseif strcmp(Mutant,'deltaPleD')%read up more on this
    param(61)=0;
elseif strcmp(Mutant,'dFtsZ')%read up more on this
    param(108)=0;
% elseif strcmp(Mutant,'deltaCtrA')
%     kdCtrA1=200;
elseif strcmp(Mutant,'deltaPodJ')
    param(37)=0;
elseif strcmp(Mutant,'deltaPdeA')
    param(64)=0;
elseif strcmp(Mutant,'divL(A601L)')
    param(96)=Pset(96)/1000; param(97)=Pset(97)*100;
elseif strcmp(Mutant,'divL(Y550F)')
%     param(96)=Pset(96)*10; param(103)=Pset(103)*10;
    param(97)=Pset(97)/10;
elseif strcmp(Mutant,'PpleC::Tn & divL(Y550F)')
    param(75)=0;
    param(121)=0.98;
%     param(96)=Pset(96)*10; param(103)=Pset(103)*10;
    param(97)=Pset(97)/10;
elseif strcmp(Mutant,'divLts')
    param(94)=param(94)/10000;
elseif strcmp(Mutant,'cckA(Y514D)')
    param(101)=Pset(101)/10000;
elseif strcmp(Mutant,'cckA(Y514D) & PdivK::Tn')
    param(101)=Pset(101)/10000;
    param(22) =Pset(22)*0.1; param(23)=Pset(23)*0.1;
% elseif strcmp(Mutant,'cdG synthetic')
%     kscdGsynthetic=0.006;
%     kscdG2=0;kscdG=0; kdcdG=0;
%     kdcdgsynthetic=0.06;
elseif strcmp(Mutant,'cckAm')
    param(10)=0;
elseif strcmp(Mutant,'EtOH_Stress')
    param(137)=1;
elseif strcmp(Mutant,'deltahdaA')
    param(114)=0;
elseif strcmp(Mutant,'WT') == 0
    error(strcat('Invalid Mutant Specification ', 'Mutant'))
end
end

