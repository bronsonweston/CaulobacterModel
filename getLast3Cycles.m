function [t1,y1,t2,y2,t3,y3,cellarrest] = getLast3Cycles(tout,yout,teout,ieout)
%[t,y,cellarrest] = getLastCycle(tout,yout,ieout)
%Takes a simulation and returns only the last cycle of the sim
t1=[];
y1=[];
t2=[];
y2=[];
t3=[];
y3=[];

for i=1:3
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
        t=tout(ind1:ind2)-teout(A(end-1));
        if teout(A(end)) < tout(end)-300 ||  t(end) > 300
            error('cell cycle arrested')
        end
        if i==1
            t1= t;
            y1= yout(ind1:ind2,:);
        elseif i==2
            t2= t;
            y2= yout(ind1:ind2,:);
        elseif i==3
            t3= t;
            y3= yout(ind1:ind2,:);
        end
        tout=tout(1:ind1+1);
        yout=yout(1:ind1+1,:);
        ieout=ieout(teout<tout(end));
        teout=teout(teout<tout(end));
        cellarrest=false;
    catch
        cellarrest=true;        
        return
    end
end

