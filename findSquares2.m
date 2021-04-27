function [score] = findSquares2(t,y, data)
yvalues=zeros(1,length(data(:,1)));
% cellcycletime= t(end);
% data(:,1)=data(:,1).*(cellcycletime/150);
for i=1:length(data(:,1))
    leftside=y(t-data(i,1)<=0);
%     rightside=y(t-data(i,1)>=0);
    yvalues(i)=leftside(end);
end

if max(data(:,2))<0.08
    d=max(data(:,2));
    n=1;
else
    d=0.25;
    n=2;
end
score = sum((abs(yvalues'-data(:,2))./d).^n./length(data(:,1)));
end

