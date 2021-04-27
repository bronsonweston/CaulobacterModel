function [] = HeatMapping(Rangex, Rangey)
% Generates heat map for Cori site d binding over a range of CtrAU and CtrA~P concentrations
% Figure 2B in main text
% Weston et al. 2021, Cell Systems
set(0,'DefaultAxesFontName', 'Arial')
titlesize=18;
ticksize=14;
labelsize=16;

%generate x and y indexes
xax = Rangex(1):Rangex(2):Rangex(3);
yax = Rangey(3):-Rangey(2):Rangey(1);
%allocate memory for figure data
Data = zeros(length(yax), length(xax));
Data2=Data;
Data3=Data;
Data4=Data;
Data5=Data;
Data6=Data;
%define parameters
kd1=0.5324*2;
kd2=0.07752/2;
kd3=0.00017/2;
%Functions to get states of DNA
DNAF=@(CtrA,CtrAP) kd1*(kd1+2*CtrA+2*CtrAP+CtrA^2/kd2+CtrAP^2/kd3+2*CtrA*CtrAP/kd2)^-1;
DNA2CtrAP2=@(CtrA,CtrAP) (DNAF(CtrA,CtrAP)*CtrAP^2)/(kd1*kd3);
DNA2CtrA2=@(CtrA,CtrAP) (DNAF(CtrA,CtrAP)*CtrA^2)/(kd1*kd2);
DNACtrA=@(CtrA,CtrAP) 2*(DNAF(CtrA,CtrAP)*CtrA)/(kd1);
DNACtrAP=@(CtrA,CtrAP) (2*DNAF(CtrA,CtrAP)*CtrAP)/(kd1);
DNACtrACtrAP=@(CtrA,CtrAP) 2*(DNAF(CtrA,CtrAP)*CtrAP*CtrA)/(kd1*kd2);

%generate data for each subplot
for i=1:length(xax)
    for j=1:length(yax) 
        ctrap=xax(i);
        ctra=yax(j);
            Data(j,i)= DNA2CtrAP2(ctra,ctrap);
            Data4(j,i)=DNACtrA(ctra,ctrap);
            Data5(j,i)=DNACtrACtrAP(ctra,ctrap);
            Data6(j,i)=DNAF(ctra,ctrap);
            Data2(j,i)=DNA2CtrA2(ctra,ctrap);
            Data3(j,i)=DNACtrAP(ctra,ctrap);
    end
    xax(i)
end
xticklabels={}; %initialize xtick labels
yticklabesl={}; %initialize ytick labels
%generate tick marks so they increase by increments of 0.1
for i=1:length(xax) 
    if rem(xax(i),max(xax)/2.2)==0
        xticklabels{i} = num2str(xax(i));
    else
        xticklabels{i} = '';
    end
end
%generate tick marks so they increase by increments of 2
for i=1:length(yax) 
    if rem(yax(i),max(yax)/2.2)==0
        yticklabels{i} = num2str(yax(i));
    else
        yticklabels{i} = '';
    end
end

%% Generate Figure
hFig=figure();
set(hFig, 'Position', [100 100 810*4/3 550])
subplot(2,4,1); 
set(gca,'fontsize',16)
colormap('jet');   % set colormap
imagesc(Data,[0 1]);        % draw image and scale colormap to values range
set(gca, 'fontsize', ticksize);
xlabel('[CtrA~P] (\muM)', 'fontsize', labelsize)
ylabel('[CtrA_U] (\muM)', 'fontsize', labelsize)
set(gca,'XTick',1:length(xax));
set(gca,'XTickLabel',xticklabels );
set(gca,'YTick',1:length(yax));
set(gca,'YTickLabel',yticklabels);
set(gca,'TickLength',[0 0])
title('[DNA:CtrA~P_2]','fontsize',titlesize)

subplot(2,4,2);
imagesc(Data2,[0 1]);        % draw image and scale colormap to values range
set(gca, 'fontsize', ticksize);
xlabel('[CtrA~P] (\muM)', 'fontsize', labelsize)
ylabel('[CtrA_U] (\muM)', 'fontsize', labelsize) 
set(gca,'XTick',1:length(xax));
set(gca,'XTickLabel',xticklabels );
set(gca,'YTick',1:length(yax));
set(gca,'YTickLabel',yticklabels);
set(gca,'TickLength',[0 0])
title('[DNA:CtrA_{U2}]', 'fontsize', titlesize)

subplot(2,4,5);
imagesc(Data3,[0 1]);        % draw image and scale colormap to values range
set(gca, 'fontsize', ticksize);
xlabel('[CtrA~P] (\muM)', 'fontsize', labelsize) 
ylabel('[CtrA_U] (\muM)', 'fontsize', labelsize) 
set(gca,'XTick',1:length(xax));
set(gca,'XTickLabel',xticklabels );
set(gca,'YTick',1:length(yax));
set(gca,'YTickLabel',yticklabels);
set(gca,'TickLength',[0 0])
title('[DNA:CtrA~P]', 'fontsize', titlesize)

subplot(2,4,6);
imagesc(Data4,[0 1]);        % draw image and scale colormap to values range
set(gca, 'fontsize', ticksize);
xlabel('[CtrA~P] (\muM)', 'fontsize', labelsize) 
ylabel('[CtrA_U] (\muM)', 'fontsize', labelsize) 
set(gca,'XTick',1:length(xax));
set(gca,'XTickLabel',xticklabels );
set(gca,'YTick',1:length(yax));
set(gca,'YTickLabel',yticklabels);
set(gca,'TickLength',[0 0])
title({'';'[DNA:CtrA_U]'}, 'fontsize', titlesize)

subplot(2,4,7);
imagesc(Data5,[0 1]);        % draw image and scale colormap to values range
set(gca, 'fontsize', ticksize);
xlabel('[CtrA~P] (\muM)', 'fontsize', labelsize) 
ylabel('[CtrA_U] (\muM)', 'fontsize', labelsize) 
set(gca,'XTick',1:length(xax));
set(gca,'XTickLabel',xticklabels );
set(gca,'YTick',1:length(yax));
set(gca,'YTickLabel',yticklabels);
set(gca,'TickLength',[0 0])
title('[DNA:CtrA~P:CtrA_U]', 'fontsize', titlesize)

subplot(2,4,3);
imagesc(Data6,[0 1]);        % draw image and scale colormap to values range
set(gca, 'fontsize', ticksize);
xlabel('[CtrA~P] (\muM)', 'fontsize', labelsize) 
ylabel('[CtrA_U] (\muM)', 'fontsize', labelsize) 
set(gca,'XTick',1:length(xax));
set(gca,'XTickLabel',xticklabels );
set(gca,'YTick',1:length(yax));
set(gca,'YTickLabel',yticklabels);
set(gca,'TickLength',[0 0])
title('[DNA]_F', 'fontsize', titlesize)

hold on
cb = colorbar(); %generates color bar
a=get(cb); %gets properties of colorbar
a.Position %gets the positon and size of the color bar
set(gca, 'fontsize', 16);
set(cb,'Position',[0.7500 0.1450 0.0300 0.7000])% To change size
print(hFig,['Figures/','HeatMap.png'],'-dpng','-r600'); %save figure 

end

