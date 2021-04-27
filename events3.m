%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This is the script for detecting specific events of the cell cycle.
% Sep. 09, 2008. Shenghua Li et al.
% This script is required by 'CauloModel.m' for running.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
function [value,isterminal,direction] = events3(t, y)
global Timer
% Locate the time when values passes through zero in a 
% increasing direction and stop integration.
%Events:
% 1) Cori hits 1 Initiation of DNA replication Swarmer (t=20, V=1.121)
% 2) Cori hits 1 Initiation of DNA replication Swarmer (t=20, V=1.121)

% 2) CcrM gene hemimethylated... [hccrM] = 1
% 3) CtrA gene hemimethylated... [hctrA] = 1
% 4) FtsZ gene hemimetnylated... [hfts] = 1 (currently not in use)
% 5) Elongation terminated... [Elong] = 0
% 6) Cell divides... [Z]=0.7 (on the way up) (t=150 for swarmer, V=
% V-(1.0103+H*0.1757)
% 7) Zring completely closed (t~132) 
% 8) PerP gene hemimeythylated

CtrA=y(1); CtrAP=y(2); DnaA=y(3); GcrA=y(4); DivK=y(5); DivKP=y(6);
I=y(7); CcrM=y(8); DivJ=y(9); CckAk=y(10); PleC=y(11); PleCk=y(12);
Elong=y(13); h_cori=y(14); h_ctra=y(15); h_ccrm=y(16); V=y(17);
ClpXP=y(18); Z=y(19); DivL=y(20); Cori=y(33);Zring=y(47);

%0.36
value=[sign(Cori - 0.99); sign(Elong - 0.25); sign(Elong - 0.36); sign(Elong - 0.615); sign(Elong-1); sign(t-(Timer+20)); sign(-1*Zring); sign(Elong - 0.75);sign(Cori - 0.999); sign(CcrM-0.65)];
isterminal=[1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
direction=[+1; +1; +1; +1; +1; +1; +1; +1; +1; +1];
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%