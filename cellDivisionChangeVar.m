function [y0] = cellDivisionChangeVar(yout,celltype, mutant, param)
y0=yout;
if strcmp(celltype, 'ST')
    celltype = 1;
elseif strcmp(celltype, 'SW')
    celltype = 0;
else
    celltype = 0;
end

K_DivL_CckA = param(100);
KclpxpCpdR = param(117);
y0(46) = 0.5;
y0(42) = 0;
y0(46) = 1;
y0(49)=y0(49)/2;
if celltype==0 %if swarmer
    if ~strcmp(mutant, 'PpleC::Tn') %CckA delocalized (and presumably DivL) in PleC mutant
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
    
else %if stalked cell
    if ~strcmp(mutant, 'PpleC::Tn')
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
    y0(25) = y0(25)*y0(17)/0.54; % [CpdR]= [CpdR]/0.54
    y0(9) = y0(9)/0.54; % [DivJ]new = [DivJ]/0.54
    y0(35)= y0(35)/0.54; % [DivJDivK]new = [DivJDivK]/0.54
    y0(36)= y0(36)/0.54; % [DivJDivK~P]new = [DivJDivK~P]/0.54
    %                 y0(17) = (0.54*y0(17)); %volume =1.08
    y0(17)=y0(17)*0.54;
end
end

