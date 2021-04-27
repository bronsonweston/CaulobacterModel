function [Pbest,BestStorage,ScoreStorage] = parameterEstimationMCMCsamples(samplesize, paramSet,paramTarget,stdv,tau,iterations1,iterations2, fname, mutants, SwarmerOn,CtrAversion)
% [Pbest,BestStorage,ScoreStorage] = parameterEstimationMCMCsamples(... 
% samplesize, paramSet,paramTarget,stdv,tau,iterations1,iterations2, ...
% saveBool, mutants, SwarmerOn)
% Returns a best parameterset

global checkmeouty0 mu
filenamelast=[];
filename=[];
% checkmeouty0 = [];
errorcount=0;
if paramTarget == 0
    paramTarget=1:length(paramSet);
end
ScoreStorage=zeros(1,iterations1*iterations2+1);
BestStorage=[];
Pbest=paramSet;
[~, ~, ~, ~, ~] = CauloSim('SW','WT', 2500,Pbest,[]);
y0=checkmeouty0;
Pcurr=Pbest;
try
    Sbest = getCost(Pbest, mutants, SwarmerOn, y0,CtrAversion)
    disp('y0=y0')
catch
    Sbest = getCost(Pbest, mutants, SwarmerOn, [],CtrAversion)
    disp('y0=[]')
    y0=[];
end
ScoreStorage(1)=Sbest;
iter=0;
iter2=0;
f=waitbar((iter+(iter2*iterations2))/(iterations1*iterations2),['Progress:', ' ', num2str(100*(iter+(iter2*iterations2))/(iterations1*iterations2)), '%']);
successPchange=0;
Plast=[];
for i=1:iterations1
    if isa(samplesize, 'char')
        paramTarg=paramTarget;
    else
        paramTarg = randsample(paramTarget,samplesize);
        
    end
    iter=0;
    while iter < iterations2
        iter=iter+1;
        if successPchange==1
            Pindex=find(Plast~=0);
            Pcurr=zeros(1,length(Plast));
            Pcurr(Pindex)=Pbest(Pindex).*(Pbest(Pindex)./Plast(Pindex));
            maxOnes = Pcurr(125:end);
            maxOnes(maxOnes>1)=1;
            Pcurr=[Pcurr(1:124),maxOnes];
            Plast=Pbest;
%             disp('success Pcurr')
        else
            Plast=Pcurr;
            Pcurr=randomPar(Pbest,paramTarg,stdv);
%             disp('standard Pcurr')
        end
        try
            tic
            [~, ~, ~, ~, ~] = CauloSim('SW','WT', 2000,Pcurr,y0.');
            y0new=checkmeouty0;
            Scurr = getCost(Pcurr, mutants,SwarmerOn, y0new,CtrAversion)
            toc
            errorcount=0;
        catch
            disp('error')
            iter=iter-1;
            errorcount=errorcount+1;
            if errorcount > 10
                [~, ~, ~, ~, ~] = CauloSim('ST','WT', 2000,Pcurr,y0.');
                y0new=checkmeouty0;
            end
            if errorcount >15
                load(filenamelast);
                Pcurr=Pbest;
                y0new=y0;
                Scurr=Sbest;
            end
            continue
        end
        Bool=MCMC(tau,Scurr,Sbest);
        if Bool==true
            successPchange=1;
            Pbest=Pcurr;
            Sbest=Scurr;
            Pbest(134)=mu;
            y0=y0new;
%             [~,~,cellarrest] = getLastCycle(tout,yout,teout,ieout);
%             if cellarrest == true
%                 y0=[];
%             end
%             BestStorage= [BestStorage,iter+1+(iter2*iterations2)];
            if ~isempty(fname) %&& Sbest < 10
                datetim=datetime('now');
                DateString = datestr(datetim);
                DateString(DateString==':')='.';
                filenamelast=filename;
                filename=[fname,'/',DateString, ' ','Svalue_', num2str(round(Sbest,4)),'.mat'];
                save(filename,'Pbest','Sbest', 'y0')
            end
        else
            successPchange=0;
        end
        waitbar((iter+(iter2*iterations2))/(iterations1*iterations2),f,['tau=',num2str(tau),', Progress',':', ' ', num2str(100*(iter+(iter2*iterations2))/(iterations1*iterations2)), '%'])
        ScoreStorage(iter+1+(iter2*iterations2))=Scurr;
    end
    iter2=iter2+1;
end
    close(f);
end

