function [ XY ] = DimerSolve( x,y,k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% a=k; b=k*y-k*x+1; c=-x;
a=1; b=-x-y-1/k; c=x.*y;
XY= (-b-(b.^2-4*a*c).^(1/2))./(2*a);
% if XY <0
%     if (x > -1e-10 && x<0) || (y > -1e-10 && y<0)
%         XY=0;
%     else
%         x
%         y
% %     disp(['X=', num2str(x), ' Y=', num2str(y), ' XY=', num2str(XY)]) 
%         error('went less than 0')
%     end
% end
end