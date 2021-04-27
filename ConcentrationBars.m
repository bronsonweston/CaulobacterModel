function  [] = ConcentrationBars(t, yset, plotvalues, labels, colors )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes hereclose all
spacings= 2:2:(2*length(plotvalues));
scale=1;
figure()
set(gcf, 'Outerposition', [100 50 1200 300]);
% colors=flipud(colors);
hold on
for j=1:length(plotvalues)
    for i=1:length(t)-1
        colorshade=[(1+(colors(j,1)-1)*yset(i,plotvalues(j))/max(yset(:,plotvalues(j))))^(1/scale),(1+(colors(j,2)-1)*yset(i,plotvalues(j))/max(yset(:,plotvalues(j))))^(1/scale),(1+(colors(j,3)-1)*yset(i,plotvalues(j))/max(yset(:,plotvalues(j))))^(1/scale)];
        colorshade(colorshade > 1)=1;
        h=fill([t(i), t(i+1), t(i+1), t(i)], ...
            [spacings(j), spacings(j), spacings(j)-1, spacings(j)-1],colorshade);
        set(h,'EdgeColor','none')
    end
end
hold off
xlim([0 t(end)])
ylim([0 spacings(end)+1])
set(gca,'XTick',0:30:t(end))
set(gca,'YTick',spacings)
set(gca,'YTickLabel',labels)
end

