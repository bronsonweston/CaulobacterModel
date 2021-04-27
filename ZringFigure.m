close all
clc
clear all

FtsA = [5 51.159;30 92.3;60 102;90 227;120 223;150 124];
FtsA(:,2)=FtsA(:,2)./max(FtsA(:,2));
FtsK= [5 22.7;30 20;60 53;90 101;120 94; 150 37];
FtsK(:,2)=FtsK(:,2)./max(FtsK(:,2));
FtsB = [5 415;30 334; 60 149;90 266; 120 442;150 468];
FtsB(:,2)=FtsB(:,2)./max(FtsB(:,2));
FtsQ = [5 101;30 79;60 54;90 241;120 378;150 157];
FtsQ(:,2)=FtsQ(:,2)./max(FtsQ(:,2));
KidO= [5 213;30 57.9;60 53;90 215.7;120 430.6;150 198.5];
KidO(:,2)=KidO(:,2)./max(KidO(:,2));
AmiC= [5 186;30 73.4;60 61.79;90 126;120 264;150 158];
AmiC(:,2)=AmiC(:,2)./max(AmiC(:,2));
GdhZ= [5 151;30 122.5;60 123;90 164;120 211;150 140];
GdhZ(:,2)=GdhZ(:,2)./max(GdhZ(:,2));

figure()
p1 = line(FtsA(:,1), FtsA(:,2),'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p1 = line(FtsK(:,1), FtsK(:,2),'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
p1 = line(FtsB(:,1), FtsB(:,2),'Color', 'b', 'LineWidth', 2, 'Linestyle', '-');
p1 = line(FtsQ(:,1), FtsQ(:,2),'Color', 'g', 'LineWidth', 2, 'Linestyle', '-');
p1 = line(KidO(:,1), KidO(:,2),'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p1 = line(AmiC(:,1), AmiC(:,2),'Color', 'c', 'LineWidth', 2, 'Linestyle', '-');
p1 = line(GdhZ(:,1), GdhZ(:,2),'Color', 'y', 'LineWidth', 2, 'Linestyle', '-');
ylabel('mRNA')
ylim([0 1.1])
xlabel('Time (min)')
set(gca,'FontSize', 14,'fontname', 'Times New Roman');
box on
legend('ftsA','ftsK','ftsB','ftsQ','kidO','amiC','gdhZ','FontAngle', 'italic')
