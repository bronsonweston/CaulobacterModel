function [CountOfIncompleteCycles] = FindCompleteCellCyclesScore(listOfMutants,Pset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
CountOfIncompleteCycles=0;
arrestListST=zeros(1,length(listOfMutants));
arrestListSW=zeros(1,length(listOfMutants));
parfor i=1:length(listOfMutants)
    mutant=string(listOfMutants(i));
    [tout, teout, ieout, yout, ~] = CauloSim('ST',mutant, 1600,Pset);
    [~,~,cellcyclearrest]=getLastCycle(tout,yout,teout,ieout);
    arrestListST(i)=cellcyclearrest;
    CountOfIncompleteCycles=CountOfIncompleteCycles+cellcyclearrest;
    [tout, teout, ieout, yout, ~] = CauloSim('SW',mutant, 1600,Pset);
    [~,~,cellcyclearrest]=getLastCycle(tout,yout,teout,ieout);
    arrestListSW(i)=cellcyclearrest;
    CountOfIncompleteCycles=CountOfIncompleteCycles+cellcyclearrest;
end

for i=1:length(listOfMutants)
    disp(strcat(string(listOfMutants(i)), ':'));
    disp(strcat('   Stalked cell cycle arrest = ', num2str(arrestListST(i))));
    disp(strcat('   Swarmer cell cycle arrest = ', num2str(arrestListSW(i))));
end
end
