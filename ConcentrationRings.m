function  [] = ConcentrationRings(t, yset, plotvalues, labels, colors, ri, Ri, events)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes hereclose all
scale=1;
figure()
set(gcf, 'Position', [1 41, 1366, 651]);
% colors=flipud(colors);
axis equal
tini=events(events(:,1)==1,2)
tend=events(events(:,1)==5,2)
tZring=events(events(:,1)==7,2)

shift= tini*2*pi/length(t)
p = linspace(pi/2+shift,-3*pi/2+shift, length(t));
p2 = linspace(pi/2+shift,-3*pi/2+shift, 10*length(t));

hold on

% yset(:,2)=yset(:,1)+yset(:,2);
for j=1:length(plotvalues)
    j
    r=ri+(Ri-ri)*(j-1);
    R=Ri+(Ri-ri)*(j-1);
    for i=1:length(t)-1
        colorshade=[(1+(colors(j,1)-1)*yset(i,plotvalues(j))/max(yset(:,plotvalues(j))))^(1/scale),(1+(colors(j,2)-1)*yset(i,plotvalues(j))/max(yset(:,plotvalues(j))))^(1/scale),(1+(colors(j,3)-1)*yset(i,plotvalues(j))/max(yset(:,plotvalues(j))))^(1/scale)];
        colorshade(colorshade > 1)=1;
        %     h=fill([p(i), p(i+1), t(i+1), t(i)], ...
        %         [spacings(j), spacings(j), spacings(j)-1, spacings(j)-1],colorshade);
        x = r*cos(p(i:i+1));
        y = r*sin(p(i:i+1));
        X = R*cos(p(i:i+1));
        Y = R*sin(p(i:i+1));
        h = patch([x flip(X)],[y flip(Y)],colorshade);
        if j==4
            if i==tini
                line([X(1) 5*(X(1)-x(1))+X(1)],[Y(1) 5*(Y(1)-y(1))+Y(1)],'color','k','linewidth',2);
            elseif i==tend
                line([X(1) 5*(X(1)-x(1))+X(1)],[Y(1) 5*(Y(1)-y(1))+Y(1)],'color','k','linewidth',2);
            elseif i==tZring
                line([X(1) 5*(X(1)-x(1))+X(1)],[Y(1) 5*(Y(1)-y(1))+Y(1)],'color','k','linewidth',2);
            end
        end
        set(h,'EdgeColor','none')
    end
    LW=1.5;
    x = r*cos(p2);
    y = r*sin(p2);
    X = R*cos(p2);
    Y = R*sin(p2);
%     L(1) = line(x,y,'color','k','linewidth',LW);
%     L(2) = line(X,Y,'color','k','linewidth',LW);
    if j==1
        L(1) = line(x,y,'color','k','linewidth',LW);
    elseif j==4
        L(2) = line(X,Y,'color','k','linewidth',LW);
    end
end

axis off
hold off
print(gcf,['Figures/','ConcentrationRings4.png'],'-dpng','-r300');
% xlim([0 t(end)])
% ylim([0 spacings(end)+1])
% set(gca,'XTick',0:30:t(end))
% set(gca,'YTick',spacings)
% set(gca,'YTickLabel',labels)
end

