function [cost] = costFun(t, y, param, mutant, mutantWeight)
%Calculates cost of simulation
if nargin <4
    mutant = 'WT';
    mutantWeight=1;
end

if strcmp(mutant, 'WT')
    cost = mutantWeight*costFunWT(t,y,param);
elseif strcmp(mutant, 'PpleC::Tn')
    cost = mutantWeight*costFunPleCko(t,y,param);
elseif strcmp(mutant, 'o')
end
end

