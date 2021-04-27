function [t,y,cellarrest] = getLastCycle(tout,yout,teout,ieout)
%[t,y,cellarrest] = getLastCycle(tout,yout,ieout)
%Takes a simulation and returns only the last cycle of the sim
try
    A = find(ieout==6);
    if A(end)-A(end-1)==1
        A=A(1:end-1);
    end
    ind1=find(tout==teout(A(end-1)));
    if isempty(ind1)&& A(end-1)==A(end-2)+1
        ind1=find(tout==teout(A(end-2)));
    end
    ind2=find(tout==teout(A(end)));
    t= tout(ind1:ind2)-teout(A(end-1));
    y= yout(ind1:ind2,:);
    if teout(A(end)) < tout(end)-300 || t(end)>300
        error('cell cycle arrested')
    end
    cellarrest=false;
catch
    t=tout;
    y=yout;
    cellarrest=true;
end

