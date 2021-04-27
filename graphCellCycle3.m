function [] = graphCellCycle3(tout, teout, ieout, yout,cellCycleIniTimes, param, celltype, TITLE)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if strcmp(celltype, 'ST')
    celltype = 1;
elseif strcmp(celltype, 'SW')
    celltype = 0;
else
    celltype = 0;
end

if nargin<8
    TITLE = [];
end

tspan = [0 tout(end)];
%% Data Points
dpCtrAP = [10*1.07 0.4; 80*1.07 0.5; 105*1.07 .8; 130*1.07 1]; %Jacobs 2003 (relative to max = 1) removed second point (t=40, P=0.15) as CtrAT suggests it is incorrect
Mcgrath = [0	12915.912; 20	9480.648; 40	2687.648; 60	515.698; 80	5769.719; 100	10582.426; 120	13602.205;140	13106.439];
Mcgrath(:,2)=Mcgrath(:,2)/max(Mcgrath(:,2));
dpCtrAT = Mcgrath;
% dpCtrAT2 = [0 0.8; 20 0; 40 0; 60 0; 80 0.5; 100 1; 120 .85; 140 0.7]; %Collier 2006 - DnaA couples DNA... (Normalized based on max expression)
dpDnaA = [0 0.78;15 1; 30 0.9; 45 0.725; 60 0.375; 75 0.325; 90 0.225; 105 0.42; 120 0.4; 135 0.6]; % Cheng and Keiler 2009 (cell cycle time 136 +- 11 min)
dpDnaA(:,1) = dpDnaA(:,1)*145/136;
dpGcrA2 = [0 0; 20 .5; 40 1; 60  0.9; 80 0.98; 100 0.58; 120 .65; 140 .75];
dpGcrA = [0 0; 20 .2; 40 .7;60  1; 80 0.8; 100 0.6; 120 .5; 140 .4]; % Holtzendorff-2004- Oscillating global regulatiors
dpPleCpole = [12*0.93 0.68; 28*0.93 0.37; 40*0.93 0.25; 60*0.93 0.4; 75*0.93 0.5; 90*0.93 0.72; 105*0.93 0.85; 120*0.93 1; 130*0.93 0.92; 150*0.93 0.72]; %Christen 2010
dpPleCtot = [0 .85/0.85; 20 0.5/0.85; 40 .28/0.85; 60 .38/0.85; 80 .45/0.85; 100 .65/0.85; 120 0.75/0.85; 140 0.8/0.85]; %Viollier 2002 - ommitted last datapoint bc it was at time point 160 and is unlikely behavior
dpPodJ1 = [0 0; 20 0; 40 .1; 60 .55; 80 1; 100 1; 120 0.42; 140 0.18]; %Chen et al.
dpPodJ2 = [0 0; 20 0.05; 40 0.1; 60 0.9; 80 1; 100 0.85; 120 0.27; 140 0.15]; %Vollier et al.
dpPerP = [0*0.93 0.6; 20*0.93 0.25; 40*0.93 0; 60*0.93 0; 95*0.93 0.5; 120*0.93 1; 130*0.93 0.8; 150*0.93 0.6]; % Shengua Li's model

dpDivJ = [0 0.364; 20 0.574; 40 0.81; 60 0.997; 80 0.924; 100 0.9; 120 1]; % Wheeler et al. 1999
dpDivJ2 = [0 0.26; 20 0.56; 40 0.86; 60 0.81; 80 0.97; 100 0.93; 120 0.86]; %Sanselicio 2015
dpDivKP = [0 0.68; 30 0.87; 60 0.87; 90 0.96; 120 1]; %Jacobs et al. 2001 (skewed from a 100 min cell cycle)
dpDivKtot = [0 0.88; 25 0.84; 50 0.79; 75 0.98; 100 0.98; 125 1]; %Jacobs et al. 2001 (skewed from a 120 min cell cycle)
%Rough estimate from Hubert Lam 2003 (0 means low or min level)
dpCckAk = [10 0.4; 40 0.3; 70 0.4; 100 1; 130 1]; %Jacobs 2003 (relative to max = 1)
dpCcrM = [5 .15; 45 0; 90 0.5; 120 1; 150 0.25]; % Grunenfelder et al. 2001 Proteomic analysis of the bacterial cell cycle (relative to max)
dpDNA = [0 0; 5 0; 10 0; 15 0.05; 20 0.18; 25 0.25; 30 0.4; 40 0.65; 50 0.75; 60 0.9; 70 0.85; 80 1; 90 0.9; 100 0.95]; % Quon et al. 1997
dpCpdRtot = [1*0.93 0.97; 25*0.93 1; 50*0.93 0.66; 75*0.93 0.497; 100*0.93 0.502; 125*0.93 0.64; 150*0.93 0.96]; %Iniesta et al. 2006 Figure 5A
dpCpdR = [1*0.93 0.42; 25*0.93 0.938; 50*0.93 0.636; 75*0.93 0.48; 100*0.93 0.07; 125*0.93 0.275; 150*0.93 0]; %Iniesta et al. 2006 Figure 5A estimated
dpcdG = [0 1/3; 20 1; 40 0.55; 60 0.35; 80 0.3; 100 0.25; 120 0.05; 140 0.35]; %Abel 2013
dpRcdA=[0 0.1; 20 0.95; 40 1; 60 0.55; 80 0.55; 100 0.5; 120 0.4; 140 0.7]; %Mcgrath 2006
dpPdeA = [0 1; 18.75 0.41; 37.5 0.158; 56.25 0.065; 75 0.072; 93.75 0.173; 112.5 0.44; 131.25 0.82; 150 0.998]; %Abel 2011 (my analysis of western)
dpFtsA = [0 250; 15 150; 30 0; 45 0; 60 355; 75 450; 90 800; 105 1700; 120 2300; 135 3200; 150 1600]; %in molecules Martin 2004
dpSciP = [0 0.87; 20 0.856; 40 0.043; 80 0; 100 0.02; 120 0.4; 140 1]; %Tan et al. 2010
%dpPleD =
%dpPleCP = [15 1; 45 0; 80 0; 135 1];
dpDNAmeth = [105 130];

for j=0:1
    cellcycletimes= teout(ieout==7);
    cellCycleEventTimes = cellCycleIniTimes;
    if j==1
        try
            A = find(ieout==6);
            if A(end)-A(end-1)==1
                A=A(1:end-1);
            end
            ind1=find(tout==teout(A(end-1)));
            ind2=find(tout==teout(A(end)));
            tout= tout(ind1:ind2)-teout(A(end-1));
            for a=length(cellCycleIniTimes):-1:1
                if cellCycleIniTimes(a)<teout(A(end))
                    cellCycleIniTimes = cellCycleIniTimes(a);
                    break
                end
            end
            cellCycleEventTimes = cellCycleIniTimes-teout(A(end-1));
            cellcycletimes= cellcycletimes-teout(A(end-1));
            tspan=[tout(1) tout(end)];
            yout= yout(ind1:ind2,:);
        catch
            continue
        end
    end
    %% calculate concentrations
    K_DivL_CckA = param(100); KclpxpCpdR = param(117); Krcdapopa=param(54);
    CckAT= yout(:,10); DivL=yout(:,34);DivLDivKP=yout(:,40); DivLDivK=yout(:,19);
    CckAcdG = yout(:,18); PopA= yout(:,27); PopAcdG2=yout(:,41);
    RcdA = yout(:,26); CpdR= yout(:,25); V= yout(:,17); 
    
    DivLT = DivL+DivLDivKP+DivLDivK;
    CckADivLT = DimerSolve(CckAT,DivLT,K_DivL_CckA);
    CckADivLDivKP = CckADivLT.*DivLDivKP./DivLT;
    CckAP = CckADivLDivKP+CckAcdG-CckADivLDivKP.*CckAcdG./CckAT;
    CckAK2 = (DivLDivK./DivLT).*CckADivLT.*(1-CckAcdG./CckAT);
    CckAK1 = CckAT-CckAK2-CckAP;
    RcdAPopA=DimerSolve(PopA+PopAcdG2, RcdA, Krcdapopa);
    ClpXP= (CpdR./(CpdR+KclpxpCpdR./V)).*RcdAPopA.*PopAcdG2./(PopAcdG2+PopA);

%     figure(j+1);
    figure();
    set(gcf, 'Name', 'Numerical Simulation of Model for Caulobacter Cell Cycle');
    set(gcf, 'Outerposition', [100 50 650 800]);
    %         set(gcf, 'Outerposition', [0 0 1280 800]);

    subplot(4, 1, 1)
        if ~isempty(TITLE)
        title(TITLE)
    end
    p1 = line(tout, yout(:, 13), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
    p2 = line(tout, yout(:, 17), 'Color', 'c', 'LineWidth', 2, 'Linestyle', '-');
    p5 = line(tout, yout(:, 33), 'Color', 'm', 'LineWidth', 2, 'Linestyle', ':');
    p1 = line(tout, yout(:, 47), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '--');
    hold on
%     if j==1
%         p4 = scatter(dpDNA(:,1)-30*celltype,dpDNA(:,2).*max(yout(:,13)), 'r');
%     end
    axis([tspan 0 2.2]);
    ax1 = gca;
    hold on
    for i=1:length(cellCycleEventTimes)
        xspot = cellCycleEventTimes(i);
        plot([xspot xspot], [0 2.5], '--k')
    end
    for i=1:length(cellcycletimes)
        xspot= cellcycletimes(i);
        plot([xspot xspot], [0 2.5], '--b')
    end
    hold off
    h = legend('Elong', 'Volume', 'Cori', 'Zring', 'Location', 'North');
    set(h, 'Interpreter', 'none', 'Box', 'on', 'Orientation', 'horizontal');
    
    subplot(4, 1, 2) %CtrA
    p1 = line(tout, yout(:, 1), 'Color', 'b', 'LineWidth', 2, 'Linestyle', ':');
    p2 = line(tout, yout(:, 2), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
    p3 = line(tout, yout(:, 1)+yout(:, 2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
    if max(yout(:,1)+yout(:,2)) > 0
        axis([tspan 0 max(yout(:,1)+yout(:,2))]);
    end
    ax1 = gca;
    hold on
    if j==1
        [~,scalarCtrAT] = findSquares(tout,yout(:, 1)+yout(:, 2), dpCtrAT);
        p4 = scatter(dpCtrAT(:,1)-30*celltype,dpCtrAT(:,2).*scalarCtrAT, 'ks');
        [~,scalarCtrAP] = findSquares(tout,yout(:, 2), dpCtrAT);
        p6 = scatter(dpCtrAP(:,1)-30*celltype,dpCtrAP(:,2)*scalarCtrAP, 'rs');
    end
    for i=1:length(cellCycleEventTimes)
        xspot = cellCycleEventTimes(i);
        plot([xspot xspot], [0 max(yout(:, 1)+yout(:, 2))*2], '--k')
    end
    for i=1:length(cellcycletimes)
        xspot= cellcycletimes(i);
        plot([xspot xspot], [0 max(yout(:, 1)+yout(:, 2))*2], '--b')
    end
    h = legend('CtrA','CtrA~P', 'CtrA Total', 'Location', 'North');
    set(h, 'Interpreter', 'none', 'Box', 'on', 'Orientation', 'horizontal');
    ylabel('C', 'Rotation', 0);
    
    hold off
    subplot(4, 1, 3)
    % p1 = line(tout, yout(:, 39), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
    p1 = line(tout, yout(:, 3), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
    p2 = line(tout, yout(:, 4), 'Color', 'g', 'LineWidth', 2, 'Linestyle', '-');
    axis([tspan 0 1.6]);
    ax1 = gca;
%     set(ax1, 'XTick', tspan(1):50:tspan(2), 'YTick', 0:0.5:1.6, 'Fontsize', 9, 'Box', 'on');
    set(h, 'Interpreter', 'none', 'Box', 'on', 'Orientation', 'horizontal');
    ylabel('D', 'Rotation', 0);
    axis([tspan 0 1.1*max([yout(:, 3);yout(:, 4)])])
    hold on
    if j==1
        [~,scalarDnaA] = findSquares(tout,yout(:, 3), dpDnaA);
        p3 = scatter(dpDnaA(:,1)-30*celltype,dpDnaA(:,2).*scalarDnaA, 'bs');
        p4 = scatter(dpGcrA2(:,1)-30*celltype,dpGcrA2(:,2).*max(yout(:,4)), 'g+');
        p5 = scatter(dpGcrA(:,1)-30*celltype,dpGcrA(:,2).*max(yout(:,4)), 'gs');
    end
    for i=1:length(cellCycleEventTimes)
        xspot = cellCycleEventTimes(i);
        plot([xspot xspot], [0 2.5], '--k')
    end
    for i=1:length(cellcycletimes)
        xspot= cellcycletimes(i);
        plot([xspot xspot], [0 2.5], '--b')
    end
    hold off
    h = legend('DnaA', 'GcrA', 'Location', 'North');

    
    subplot(4, 1, 4)
    p1 = line(tout, yout(:, 14), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
    p2 = line(tout, yout(:, 15), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
    p3 = line(tout, yout(:, 16), 'Color', 'g', 'LineWidth', 2, 'Linestyle', ':');
    p6 = line(tout, yout(:, 24), 'Color', 'c', 'LineWidth', 2, 'Linestyle', '-');
    p6 = line(tout, yout(:, 45), 'Color', 'y', 'LineWidth', 2, 'Linestyle', '-');
    p5 = line(tout, yout(:, 8), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-.');
    axis([tspan 0 1.6]);
    ax1 = gca;
    set(ax1, 'XTick', tspan(1):50:tspan(2), 'YTick', 0:0.5:1.6, 'Fontsize', 9, 'Box', 'on');
     ylabel('B', 'Rotation', 0);
    hold on
    if j==1
        p6 = line([dpDNAmeth(1)-30*celltype dpDNAmeth(1)-30*celltype], [0 2], 'Color', 'm', 'LineWidth', 2, 'Linestyle', '--');
        p7 = line([dpDNAmeth(2)-30*celltype dpDNAmeth(2)-30*celltype], [0 2], 'Color', 'm', 'LineWidth', 2, 'Linestyle', '--');
        p8 = scatter(dpCcrM(:,1)-30*celltype,dpCcrM(:,2).*max(yout(:,8)), 'm');
    end
    for i=1:length(cellCycleEventTimes)
        xspot = cellCycleEventTimes(i);
        plot([xspot xspot], [0 2.5], '--k')
    end
    for i=1:length(cellcycletimes)
        xspot= cellcycletimes(i);
        plot([xspot xspot], [0 2.5], '--b')
    end
    h = legend('hCori', 'hctrA', 'hccrM', 'hperP', 'hfts', 'CcrM', 'Location', 'North');
    set(h, 'Interpreter', 'none', 'Box', 'on', 'Orientation', 'horizontal');
    

end
end

