function [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype, mutant, TimeOfRun,Pset, y0)
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% S[tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype, mutant, TimeOfRun,param)
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% celltype determines which cell cycle (swarmer or stalked) we are tracking.
if length(Pset) < 139
    Pset=[Pset,1,1,1,0.98];
end
Pset(108)=Pset(108)*Pset(32)/Pset(33);
Pset(32)=1;
Pset(33)=1;
    

param=Pset;
global  Timer checkmeouty0 mu
K_DivL_CckA = param(100);
KclpxpCpdR = param(117);

if nargin < 1
    celltype = 0;
    TimeOfRun = 1500;
    mutant = 'WT';
    Y0=y0;
end

if strcmp(celltype, 'ST')
    celltype = 1;
elseif strcmp(celltype, 'SW')
    celltype = 0;
else
    celltype = 0;
end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%
%% Initial value of model variables, for a newborn, wild-type swarner cell
if nargin < 5 || isempty(y0)
    % if isempty(y0)
    y0 = zeros(49, 1);
    
    y0(1) = 0.01; % [CtrA]
    y0(2) = 0.03;  %[CtrA~P]
    y0(3) = 0.3; % [DnaAtot]
    y0(48) = 0.28; % [DnaAatp]
    y0(4) = 0.18; %[GcrA]
    y0(5) = 0.05;  %[DivK]
    y0(6) = 0.2;    %[DivK~P]
    y0(7) = 0;  % SciP
    y0(8) = 0.05; %[CcrM]
    y0(9) = 0.004;    %DivJ
    y0(10) = 0.025;   %CckAT
    y0(11) = 0.004;  %PleC
    y0(12) = 0.048;    %PleCk
    y0(13) = 0;    %[Elong](Elongation)
    y0(14) = 1;    %[hCori]
    y0(15) = 1;    %[hctrA]
    y0(16) = 1;    %[hccrM]
    y0(17) = 1.01;  %Volume (V)
    y0(18) = 0.018;    %CckA:cdG
    y0(19) = 0.03;   % DivLDivK
    y0(20) = 1;    %X
    y0(21) = 0.043; %PleCpole
    y0(22) = 0.75; %PodJ;
    y0(23) = 0.156; %PerP
    y0(24) = 1;  %[hperp] (hemimethylated PerP)
    y0(25) = 0.4; %CpdR
    y0(26) = 0.03; %RcdA
    y0(27) = 0.007; %PopA
    y0(28)=0; %PdeA
    y0(29)=0.00055; %PleD
    y0(30)=0.002; %PleD~P
    y0(31)=0.11; %cdG
    y0(32)=0; %RepSwitch
    y0(33)=0; %Cori
    y0(34)=0.04; %DivL
    y0(35)=0.003; %DivJDivK
    y0(36)=0; %DivJDivKP
    y0(37)=0.038; %CpdRP
    y0(38)=0.0001; %PleDcdG2
    y0(39)=0.0004; %PleDPcdG2
    y0(40)=0.08; %DivLDivKP
    y0(41)=0.06; %PopAcdG2
    y0(42)= 0; %MipZswitch=y(42);
    y0(43)= 0; %DgcBcdG;
    y0(44)= 0.17; %FtsA= y(44);
    y0(45)= 1; %h_plec = y(45);
    y0(46)=0; %FtsAdeg
    y0(47)=1; %Zring
    y0(49)= 1; %ChromosomeCount
    Timer=-21;
    % end of initial values
    % y0=[0.0850433016925143,0.0502590200033401,0.360176399312856,0.214672809281808,0.0136343522382722,0.135266173060524,0.0332756339465719,0.0417892748762724,0.00606274238180525,0.0205677168092387,0.00163739733070039,0.0311385556565151,0,1,4.81997286409898e-05,4.81997286409898e-05,1.12039344494061,0.0136473497565798,0.00986550228694938,1.89744820746202,0.0125732826518279,0.00103510216830892,0.176441531810805,4.81997286409898e-05,0.327644197875673,0.0275227451576489,0.00399293859304769,3.03688193163146e-05,0.000166772315253276,0.00173052903061077,0.0980463379904026,1,0,0.0301754953604212,0.00133887509839787,0.000151014031652913,0.0732852281477257,3.87223945079591e-05,0.000401806686121587,0.0875185193849886,0.0417869339919082,0,0.00148957006024837,0.150174931355794,4.81997286409898e-05,4.67931709242343e-47,1].';
else
    Timer=-21;
end
Size=size(y0);
if Size(1)<Size(2)
    y0=y0.';
end

Y0=y0;
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Integration parameters
tspan = [0 TimeOfRun];
tstart = tspan(1);
tfinal = tspan(2);
% [T, Y] = ode15s(@Caulobacter2, [tstart tfinal], y0);


% end of integration parameters
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Numerical Integration

tout = tstart;
yout = y0.';
teout = [];
yeout = [];
ieout = [];
storageOn =0;
% if ~isempty(y0storage) && celltype == 1
%     y0=y0storage;
%     yout = y0;
%     storageOn=0;
% end

% loop for continous simulation of the differential equations.
% i is used to counted the interrupt points through the simulation. It can be
% changed to larger to satisfy longer time simulation once more interrupt
% points are detected.
timeofcycle=[];
tinitcyc=tstart;
cellCycleEventTimes=[];
count=1;
cellCycleIniTimes = [];
warning('error', 'MATLAB:ode15s:IntegrationTolNotMet')
Mutant = 'WT';
IE=0;
% if isempty(mu)
%     MU=0.0054;
% else
%     Pset(134)=mu;
% end
tol=3.0e-4;
% options = odeset('Events', @events3, 'RelTol', tol, 'AbsTol', tol/10);
% options2 = odeset('Events', @events3, 'RelTol', tol/100, 'AbsTol', tol/1000);
options = odeset('Events', @events3, 'RelTol', 3.0e-8, 'AbsTol', 3.0e-8);
options2 = odeset('Events', @events5, 'RelTol', 3e-9, 'AbsTol', 3e-10);
% options3 = odeset('Events', @events3, 'RelTol', 3e-10, 'AbsTol', 3e-11);
MU=Pset(134);
while tstart < tfinal
    %     tstart
    %     y0
    tfinal2=tfinal;
    %     try
    %     if tstart < 700
    %         y0
    %         try
    %             Difference = y0-Y(end,:)
    %         catch
    %         end
    %         tstar
    %     if ~isempty(timeofcycle) && IE==1
    %         Mutant=mutant;
    %     end
    if ~isempty(timeofcycle) && IE==7
        Mutant=mutant;
        MU=log(2/y0(17))/(timeofcycle(end));
        if MU < 0.0038
            MU=0.0038;
        elseif MU > 0.007
            MU=0.007;
        end
        Pset(134)=MU;
        mu=MU;
        param = getMutantParam(Mutant,Pset);
    end
%         y0(abs(y0)==0)=1e-10;
    try
        [T, Y, TE, YE, IE] = ode15s(@(t,y) cauloODE(t,y,param, Mutant), [tstart tfinal], y0, options);
    catch
        try
            [T, Y, TE, YE, IE] = ode15s(@(t,y) cauloODE(t,y,param, Mutant), [tstart tfinal], y0, options2);
        catch
            [tout, teout, ieout, yout, cellCycleIniTimes] = CauloSimFineStep(celltype, mutant, TimeOfRun,Pset, Y0);
            return
        end
    end
    %         Y(end,:)
    %     else
    %         [T, Y, TE, YE, IE] = ode15s(@Caulobacter2DecreasedNutrients, [tstart tfinal], y0, options);
    %     end
    %     catch
    %         try
    %             [T, Y, TE, YE, IE] = ode15s(@Caulobacter2temp, [tstart tstart+1], y0, options);
    %         catch
    %             'ERROR AGHHHHH'
    %             tstart
    %             y0
    %             break
    %             for i=1:length(y0)
    %                 if y0(i)~=0
    %                     y0(i)= y0(i)+0.1;
    %                 end
    %             end
    %         end
    %     end
    % Accumulate output.
    nt = length(T);
    tout = [tout;T(2:nt)];
%     Y
    yout = [yout;Y(2:nt, :)];
    teout = [teout;TE];
    yeout = [yeout;YE];
    ieout = [ieout;IE];
    
    % Refresh the initial conditions for differential equations once the interrupt
    % point is detected.
    y0 = Y(nt, :);
    
    if isscalar(IE) == 0
        IE = 0;
    end
    
    switch (IE)
        
        case 1 % 1) Initiation of DNA replication
            %             if celltype==0 || count==1
            if y0(32)==1 && y0(13)<0.1
            else
                y0(49)=y0(49)*2;
            end
            if y0(8) < 0.65
                y0(14) = 0.5;
            end
            y0(46)=0;
            y0(32) = 1;
            y0(33)=0;
            cellCycleIniTimes=[cellCycleIniTimes,T(nt)];
            %             checkmeouty0=y0;
            %             end
        case 2
            if y0(8) < 0.65
                y0(16) = 0.5;
            end
        case 3
            if y0(8) < 0.65
                y0(15) = 0.5;
            end
        case 4
            if y0(8) < 0.65
                y0(45) = 0.5;
            end
        case 5 %Elongation terminated
            y0(13) = 0;
            y0(32) = 0;
            y0(42) = 1;
            %             if strcmp(mutant,'WT') && celltype==0
            %                 checkmeouty0=y0;
            %             end
        case 6 %cell division
            count=count+1;
            timeofcycle=[timeofcycle;T(nt)-tinitcyc];
            tinitcyc=T(nt);
            %             cellCycleEventTimes=[cellCycleEventTimes,T(nt)];
            y0(47)=1;
            if strcmp(mutant,'WT')
                checkmeouty0=y0;
            end
%             y0(46)=0;
        case 7 % Zring completely closed
            y0(42) = 0;
            y0(46) = 0.5;
            y0(49)=y0(49)/2;
            Timer = T(nt);
            %             if storageOn == 1
            %                 y0storage= Y(end-2, :);
            %             end
            if celltype==0 %if swarmer
                if y0(21)>0.01
                    %CckA delocalized (and presumably DivL) in PleC mutant
                    CckADivLT = DimerSolve(y0(10),y0(34)+y0(40)+y0(19),K_DivL_CckA);
                    CckAf=y0(10)-CckADivLT;
                    DivLDivKP=y0(40);
                    DivLT = y0(34)+y0(40)+y0(19);
                    CckADivLDivKP = CckADivLT*DivLDivKP/DivLT;
                    y0(18)=(y0(18)/y0(10))*(CckAf+(y0(10)-CckAf-CckADivLDivKP)*y0(17)/(0.46*y0(17))); % [CckA:cdG]= ([CckA:cdG]/[CckA])*([CckA]-[CckA:DivLT]+([CckADivLT]-CckADivLDivKP)))
                    y0(10) = (y0(10)-CckAf-CckADivLDivKP)*y0(17)/(0.46*y0(17))+CckAf; %CCkAT
                    y0(34)=y0(34)*y0(17)/(0.46*y0(17)); %DivL (All of it at swarmer pole)
                    y0(40)=0; %DivLDivK~P (All of it is at stalked pole)
                    y0(19)=y0(19)*y0(17)/(0.46*y0(17));%DivLDivK (All of it at swarmer pole)
                end
                a = (y0(11)*y0(21)/(y0(11)+y0(12)))*y0(17)/(0.46*y0(17)) + ((y0(11)+ y0(12)-y0(21))*y0(11)/(y0(11)+y0(12)));  %PleC
                b = (y0(12)*y0(21)/(y0(11)+y0(12)))*y0(17)/(0.46*y0(17)) + ((y0(11)+ y0(12)-y0(21))*y0(12)/(y0(11)+y0(12)));    %PleCk
                y0(11)=a; % [PleC]new = [PleC]*[PleCpole]/(0.46*[PleCtot]) +  ([PleCtot]-[PleCpole])*[PleC]/[PleCtot]
                y0(12)=b; % [PleCk]new = [PleCk]*[PleCpole]/(0.46*[PleCtot]) +  ([PleCtot]-[PleCpole])*[PleCk]/[PleCtot]
                y0(21) = y0(21)*y0(17)/(0.46*y0(17)); % [PleCpole]new = [PleCpole]/0.46
                y0(26) = y0(26)-((y0(25)/(y0(25)+KclpxpCpdR/y0(17)))*y0(26)); % RcdA=RcdA-RcdA*(CpdR/(CpdR+KclpXPcpdR)) %RcdA remove fraction bound to CpdR at stalk pole
                y0(30)=(1/10)*y0(30); % [PleD~P]new= (1/10)*[PleD~P]  %PleD~P (9/10th is localized at stalked pole)
                y0(39)= (1/10)*y0(39); % [PleD~P:cdG2]= (1/10)*[PleD~P:cdG2] %PleDPcdG2 (9/10th is localized at stalked pole)
                y0(9) = 0; % [DivJ]=0
                y0(35)= 0; % [DivJDivK]=0
                y0(36)= 0; % [DivJDivK~P]=0
                y0(20)=0.01; %X
                y0(25) = 0; % [CpdR]=0
                y0(17)=y0(17)*0.46; % [V]new= V*0.46
                %                 y0(17)=0.92;
            else %if stalked cell
                %                 CckADivL = DimerSolve(y0(10),y0(34)+y0(40)+y0(19),K_DivL_CckA);
                %                 CckAf=y0(10)-CckADivL;
                %                 y0(18)=(y0(18)/y0(10))*(CckAf+0.3*CckADivL*y0(17)/(0.54*y0(17)));
                %                 y0(10) = CckAf+0.3*CckADivL*y0(17)/(0.54*y0(17)); %CCkAT
                %                 y0(34)=y0(34)*0.3*y0(17)/(0.54*y0(17)); %DivL (30% goes to pole as an arbitrary estimate)
                %                 y0(40)=y0(40)*0.3*y0(17)/(0.54*y0(17)); %DivLDivK~P (30% goes to pole as an arbitrary estimate)
                %                 y0(19)=y0(19)*0.3*y0(17)/(0.54*y0(17)); %DivLDivK~P (30% goes to pole as an arbitrary estimate)
                if y0(21)>0.01
                    CckADivLT = DimerSolve(y0(10),y0(34)+y0(40)+y0(19),K_DivL_CckA);
                    CckAf=y0(10)-CckADivLT;
                    DivLDivKP=y0(40);
                    DivLT = y0(34)+y0(40)+y0(19);
                    CckADivLDivKP = CckADivLT*DivLDivKP/DivLT;
                    y0(18)=(y0(18)/y0(10))*(CckAf+(CckADivLDivKP)*y0(17)/(0.54*y0(17)));
                    y0(10) = (CckADivLDivKP)*y0(17)/(0.54*y0(17))+CckAf; %CCkAT
                    y0(34)=0; %DivL (All of it at swarmer pole)
                    y0(40)=y0(40)*y0(17)/(0.54*y0(17)); %DivLDivK~P (All of it is at stalked pole)
                    y0(19)=0;%DivLDivK (All of it at swarmer pole)
                end
                a = ((y0(11)+ y0(12)-y0(21))*y0(11)/(y0(11)+y0(12)));  %PleC = 0
                b = ((y0(11)+ y0(12)-y0(21))*y0(12)/(y0(11)+y0(012)));  %PleCk = 0
                y0(11) = a; % [PleC]new = ([PleCtot]-[PleCpole])*[PleC]/[PleCtot]
                y0(12) = b; % [PleCk]new = ([PleCtot]-[PleCpole])*[PleCk]/[PleCtot]
                y0(30)=(9/10)*y0(30)*y0(17)/(0.54*y0(17))+(1/10)*y0(30); % [PleD~P]= (9/10)*[PleD~P]/0.54+(1/10)*[PleD~P]      %PleD~P (9/10th is localized)
                y0(39)= (9/10)*y0(39)*y0(17)/(0.54*y0(17))+(1/10)*y0(39); % [PleDPcdG2] = (9/10)*[PleD~P]/0.54+(1/10)*[PleD~P]     %PleDPcdG2 (9/10th is localized)
                y0(26) = y0(26)-(y0(25)/(y0(25)+KclpxpCpdR/y0(17)))*y0(26)+(y0(25)/(y0(25)+KclpxpCpdR/y0(17)))*y0(26)/0.54;   % [RcdA]new=[RcdA]+[RcdA]*(CpdR/(0.54*(CpdR+KclpXPcpdR))  %RcdA add fraction bound to CpdR at stalk pole
                y0(21) = 0; %PleCpole
                y0(25) = y0(25)/0.54; % [CpdR]= [CpdR]/0.54
                y0(9) = y0(9)/0.54; % [DivJ]new = [DivJ]/0.54
                y0(35)= y0(35)/0.54; % [DivJDivK]new = [DivJDivK]/0.54
                y0(36)= y0(36)/0.54; % [DivJDivK~P]new = [DivJDivK~P]/0.54
                %                 y0(17) = (0.54*y0(17)); %volume =1.08
                y0(17)=y0(17)*0.54;
                
            end
%                         checkmeouty0 = y0;
        case 8 %PerP hemimethylated
            if y0(8) < 0.65
                y0(24) = 0.5;
            end
        case 9 %  Initiation of DNA replication in swarmer cell BACKUP
            if y0(32)==1 && y0(13)<0.1
            else
                y0(49)=y0(49)*2;
            end
            y0(32) = 1;
            if y0(8) < 0.65
                y0(14) = 0.5;
            end
            y0(33)=0;
            y0(46)= 0;
            cellCycleIniTimes=[cellCycleIniTimes,T(nt)];
        case 10 %  methylation
            y0(14)=1;
            y0(15)=1;
            y0(16)=1;
            y0(24)=1;
            y0(45)=1;
    end
    
    tstart = T(nt);
    if tstart >= tfinal
        break;
    end
    
end
stuff = [teout,ieout];
% timeofcycle

%end of integration
% [tout, yout, teout, yeout, ieout ] = CauloSim(y0storage.', 0, 500, 'stalked');
% yout
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%% Save the data
%save('swstmodelwtsim30.txt', 'tout', 'yout', '-ASCII');
%Save the initial data at t  = 270 min
%y_initial = interp1(tout, yout, 442);
%save('swstmodelwtsim445.txt', 'y_initial', '-ASCII');

end
