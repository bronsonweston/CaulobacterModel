function [ XY ] = DimerSolve2( x,y,k )

a=1; b=-x-y-1/k; c=x*y;
XY= (-b-(b^2-4*a*c)^(1/2))/(2*a);
end