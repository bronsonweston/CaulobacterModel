function [Pbest,BestStorage,ScoreStorage] = parameterEstimationMCMC(paramSet,paramTarget,stdv,beta,iterations,saveBool, mutants)
%[Pbest,BestStorage,ScoreStorage] = parameterEstimationMCMC(paramSet,paramTarget,stdv,beta,iterations)
% Returns a best parameterset
if paramTarget == 0
    paramTarget=1:length(paramSet);
end
ScoreStorage=zeros(1,iterations+1);
BestStorage=[1];
Pbest=paramSet;
Sbest = getCost(Pbest, mutants);
ScoreStorage(1)=Sbest;
iter=0;
f=waitbar(iter/iterations,['Progress:', ' ', num2str(100*iter/iterations), '%']);

while iter < iterations
    iter=iter+1;
    Pcurr=randomPar(Pbest,paramTarget,stdv);
    try
        Scurr = getCost(Pcurr, mutants)
    catch
        iter=iter-1;
        continue
    end
    Bool=MCMC(beta,Scurr,Sbest);
    if Bool==true
        Pbest=Pcurr;
        Sbest=Scurr;
        BestStorage= [BestStorage,iter+1];
        if saveBool == true
            datetim=datetime('now');
            DateString = datestr(datetim);
            DateString(DateString==':')='.';
            filename=['ParameterSets\',DateString, ' ','Pvalue_', num2str(round(Sbest,4)),'.mat'];
            save(filename,'Pbest','Sbest')
        end
    end
    waitbar(iter/iterations,f,['Progress:', ' ', num2str(100*iter/iterations), '%'])
    ScoreStorage(iter+1)=Scurr;
end
close(f);
end

