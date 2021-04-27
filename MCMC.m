function [Bool] = MCMC(tau,cost,mincost)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if rand() < exp((mincost-cost)/tau)
    Bool= true;
else
    Bool= false;
end     