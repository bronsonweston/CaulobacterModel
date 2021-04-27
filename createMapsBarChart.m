function [] = createMapsBarChart()
% generates bar charts corresponding to figure 5A,5B and 5E of main text.
% Weston et al. 2021, Cell Systems
close all %close any open figures
allMutants = {'WT','CtrAdelta3omega','CtrAD51E','CtrAD51Edelta3omega','SM921','deltaGcrA','deltaccrM', ...
    'LS1','PpleC::Tn','divLts','divL(A601L)','divL(Y550F)', 'PpleC::Tn & divL(Y550F)', ...
    'deltadivJ','PpleC::Tn & deltadivJ','deltaSpmX','PpleC::Tn & deltaSpmX', 'PdivK::Tn', ...
    'divKcs','divKxyl','cdG0','PdivK::Tn & cdG0','deltapopA','deltaPleD','deltapopA & PdivK::Tn' ...
    'deltapopA & deltaPleD','deltadgcB','deltaPdeA','deltapopA & deltaPleD & PdivK::Tn', ...
    'cckA(Y514D)','cckA(Y514D) & PdivK::Tn','deltaPodJ','deltahdaA','deltaPleD & PdivK::Tn', 'deltaRcdA'}; %define list of all mutants
ArrestList=[0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0]; %define experimental observations of cell cycle arrest for each entry in "allMutants". 0 indicates viable, 1 indicates arrest
mutantArrest=allMutants(ArrestList==1); %gets indexes of arrested strains
mutantViable=allMutants(ArrestList==0); %gets indexes of viable strains
arrestind=[]; 
viableind=[];
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
set(0,'DefaultLegendAutoUpdate','off')
Maps=load('ParamCatalog\Maps.mat'); %load simulation data
mutantlist={'WT','deltaPdeA','PpleC::Tn', 'deltapopA & PdivK::Tn','deltaccrM','deltaGcrA','deltaPleD','divL(A601L)','LS1','CtrAdelta3omega','cdG0','CtrAD51E'}; %define list of mutants used in parameterization algorithm (for color coding purposes only)
indexes=fieldnames(Maps); %get indexes of each type of simulation in Maps
xlabels = keys(Maps.(indexes{1})); %get list of all strains in Maps
idx = strcmp(xlabels,'ParamChange');
xlabels(idx) = []; %remove 'ParamChange' item from xlabels
yvalues=zeros(numel(xlabels),numel(indexes)); %allocate memory space for failure rate data
for map=1:numel(indexes) %incorporate failure data into yvalues
    currentMap=Maps.(indexes{map});
    for key=1:numel(xlabels)
        yvalues(key,map)=currentMap(char(xlabels(key)))/currentMap('ParameterCount');
    end
end
yvalues(strcmp(xlabels,'ParameterCount'),:)=[]; %remove entry on 'ParameterCount'
xlabels(strcmp(xlabels,'ParameterCount'))=[];  %remove entry on 'ParameterCount'
 
%% This section modifies mutant labels to look professional for figures
% This section is not commented in detail, as the contents are not deemed 
% important for understanding the functions
longestWord=0;
for i=1:numel(xlabels)
    ct=[];
    mut=char(xlabels(i));
    strings=strsplit(mut,'_');
    try
        ct = char(strings(2));
    catch
    end
    mut=char(strings(1));
    if any(strcmp(mutantArrest,mut))
        arrestind=[arrestind,i] %notes indexes within xlabel of strains that were observed to be arrested in experiments
    elseif any(strcmp(mutantViable,mut))
        viableind=[viableind,i] %notes indexes within xlabel of strains that were observed to be viable in experiments
    else
        disp(mut)
    end
    
    
    if any(strcmp(mut,mutantlist))
        mut=['\color{red}',mut];
    end
    
    strings=strsplit(mut,'delta');

    if length(strings)>1
        mut = char(strings(1));
        for j=2:numel(strings)
            word=char(strings(j));
            word(1)=lower(word(1));
            mut=[mut,'\Delta',word];
        end
    end
    strings=strsplit(mut,'omega');
    if length(strings)>1
        mut = char(strings(1));
        for j=2:numel(strings)
            word=char(strings(j));
            mut=[mut,'\Omega',word];
        end
    end
    
    strings=strsplit(mut,'divKxyl');
    if length(strings)>1
        mut = char(strings(1));
        for j=2:numel(strings)
            word=char(strings(j));
            mut=[mut,'P{\itxyl}-divK',word];
        end
    end
    
    Names={'CtrA','gcrA','ccrM','pleC','divL','divJ','spmX','divK','popA','pleD','dgcB','pdeA','cckA','podJ','hdaA','rcdA','3','A601L','Y514D','Y550F','D51E'};
    for n=1:length(Names)
        name=Names{n};
        strings=strsplit(mut,name);
        if n==1
            name='ctrA';
        end
        if length(strings)>1
            mut = char(strings(1));
            for j=2:numel(strings)
                word=char(strings(j));
                mut=[mut,'{\it',name,'}',word];
            end
        end
    end
    
    xlabels{i}=[mut, ' ',ct];
    if length(xlabels{i})>longestWord
        longestWord=length(xlabels{i});
    end
end
yoriginal=yvalues;

for i=1:numel(xlabels)
    word=xlabels{i};
    if contains(word,'A601L') %move data for strain with leter code 'A601L'
        xlabelA=xlabels(i); xlabelB=xlabels(i+1);
        yvalueA=yvalues(i,:); yvalueB=yvalues(i+1,:);
        xlabelC=xlabels(i+2); xlabelD=xlabels(i+3);
        yvalueC=yvalues(i+2,:); yvalueD=yvalues(i+3,:);
        xlabels(i)=xlabelC; xlabels(i+1)=xlabelD;
        xlabels(i+2)=xlabelA; xlabels(i+3)=xlabelB;
        yvalues(i,:)=yvalueC; yvalues(i+1,:)=yvalueD;
        yvalues(i+2,:)=yvalueA; yvalues(i+3,:)=yvalueB;
        break
    end
end
for i=1:numel(xlabels)
    word=xlabels{i};
    if contains(word,'xyl') %move data for strain with leter code 'xyl'
        xlabelA=xlabels(i); xlabelB=xlabels(i+1);
        yvalueA=yvalues(i,:); yvalueB=yvalues(i+1,:);
        xlabelC=xlabels(i+2); xlabelD=xlabels(i+3);
        yvalueC=yvalues(i+2,:); yvalueD=yvalues(i+3,:);
        xlabels(i)=xlabelC; xlabels(i+1)=xlabelD;
        xlabels(i+2)=xlabelA; xlabels(i+3)=xlabelB;
        yvalues(i,:)=yvalueC; yvalues(i+1,:)=yvalueD;
        yvalues(i+2,:)=yvalueA; yvalues(i+3,:)=yvalueB;
        break
    end
end
for i=1:numel(xlabels)
    word=xlabels{i};
    if longestWord==length(word)
        xlabelA=xlabels(i); xlabelB=xlabels(i+1);
        yvalueA=yvalues(i,:); yvalueB=yvalues(i+1,:);
        xlabels([i,i+1])=[]; yvalues([i,i+1],:)=[];
        pos=length(xlabels)/2;
        xlabels=[xlabels(1:pos),xlabelA,xlabelB,xlabels(pos+1:end)];
        yvalues=[yvalues(1:pos,:);yvalueA;yvalueB;yvalues(pos+1:end,:)];
        break
    end
end
%% This section produces figures and is commented

%Create Copies of yvalues and xlabels
yvalues2=yvalues; 
xlabels2=xlabels;
%Delete all entries of yvalues and xlabels that fail less than 25% of the time
for i=length(xlabels2):-1:1
    if all(yvalues2(i,:)<=0.25)
        yvalues(i,:)=[];
        xlabels(i)=[];
    end
end

% Generate figure on F/T fraction part 1 (not used in publication, but same information as Figure 5A)
fig=figure();
set(gcf,'position',[186 142 456 513]);
marka=length(xlabels2)/2;
xlabels2a=xlabels2(1:marka);
xlabels2b=xlabels2(marka+1:end);
bg= barh(1:length(xlabels2a),yvalues2(1:length(xlabels2a),[3,1,2]));
bg(1).FaceColor='y';
bg(1).EdgeColor='y';
bg(3).FaceColor='b';
bg(3).EdgeColor='b';
bg(2).FaceColor='g';
bg(2).EdgeColor='g';
str=xlabels2a(1:length(xlabels2)/2);
set(gca, 'YtickLabel',str, 'Ytick',1:numel(str))
ylim=get(gca,'ylim');
hold on
plot([0.25 0.25], ylim,':k','linewidth',1.5);
xlabel('F/T Fraction')
h=gca; h.YAxis.TickLength = [0 0];

% Generate figure on F/T fraction part 2 (not used in publication, but same information as Figure 5A)
fig=figure();
set(gcf,'position',[186 142 456 513]);
bg=barh(1:length(xlabels2b),yvalues2(length(xlabels2a)+1:end,[3,1,2]));
bg(1).FaceColor='y';
bg(1).EdgeColor='y';
bg(3).FaceColor='b';
bg(3).EdgeColor='b';
bg(2).FaceColor='g';
bg(2).EdgeColor='g';
str=xlabels2b;
set(gca, 'YtickLabel',str, 'Ytick',1:numel(str))
ylim=get(gca,'ylim');
hold on
ax=gca;
legend(ax.Children.',{'Quick','Slow','Cori^-'},'Position',[0.7190 0.8160 0.1798 0.1033]);
plot([0.25 0.25], ylim,':k','linewidth',1.5);
xlabel('F/T Fraction')
h=gca; h.YAxis.TickLength = [0 0];

% Generate figure on F/T fraction, does not produce Matlab figure but saves as png (Figure 5A)
fig=figure();
set(fig, 'visible', 'off')
set(gcf, 'PaperPosition', [0 0 5 11])    % can be bigger than screen 
set(gcf, 'PaperSize', [11 11])    % Same, but for PDF output
bg= barh(1:length(xlabels2),yvalues2(1:length(xlabels2),[3,1,2]));
bg(1).FaceColor='y';
bg(1).EdgeColor='y';
bg(3).FaceColor='b';
bg(3).EdgeColor='b';
bg(2).FaceColor='g';
bg(2).EdgeColor='g';
str=xlabels2(1:length(xlabels2));
set(gca, 'YtickLabel',str, 'Ytick',1:numel(str))
ylim=get(gca,'ylim');
ax=gca;
legstr={['\fontsize{8}Q','\fontsize{6}UICK'],['\fontsize{8}S','\fontsize{6}LOW'],['\fontsize{8}C','\fontsize{6}ORI','\fontsize{8}^-']};
legend(ax.Children.',legstr);
leg=get(gca,'legend')
pos=get(leg,'position');
set(leg,'position',[pos(1)/1.05 pos(2)*1.05 pos(3) pos(4)])
hold on
plot([0.25 0.25], ylim,':k','linewidth',1.5);
xlabel('F/T Fraction')
h=gca; h.YAxis.TickLength = [0 0];
print(fig,['Figures/','FoverTfraction.png'],'-dpng','-r600');  

% Generate figure on Average strain viability (Figure 5E)
fig=figure()
falerate=0.25;
set(gcf,'position',[   100   284   498*1.1   243])
viable=zeros(1,length(indexes));
for i=1:length(indexes) %get fraction of viable simulations
    v=yoriginal(:,i);%v equals failure rate data for parameter collection 
    v(viableind)=1-yoriginal(viableind,i); %failures on viable strains are inviable. Thus, viable count is the inverse for these strains
    viable(i)=sum(v)/length(v); %viability fraction for paremeter collection i
end
a=100*(length(ArrestList)-sum(ArrestList))/length(ArrestList);
bar(1:1.5:10,[a,100*(viable([1,4,2,6,3,5]))],'LineWidth',1)
set(gca, 'ylim', [20 80])
str={'\fontsize{10}Exp \newlineObs',['\fontsize{10}S','\fontsize{8.5}LOW'],['      ','\fontsize{10}S','\fontsize{8.5}LOW', '\newline', '\fontsize{10}-CtrA_U:{\itCori}'],...
    ['\fontsize{10}Q','\fontsize{8.5}UICK'], ['    ','\fontsize{10}Q','\fontsize{8.5}UICK', '\newline', '\fontsize{10}-CtrA_U:{\itCori}'], ...
    ['\fontsize{10}C','\fontsize{8.5}ORI','\fontsize{10}^-'],['    ','\fontsize{10}C','\fontsize{8.5}ORI','\fontsize{10}^-', '\newline', '\fontsize{10}+CtrA_U:{\itCori}']};
set(gca, 'XtickLabel',str, 'Xtick',1:1.5:(1.5*(numel(str)-1)+1))
ylabel('Avg. Strain Viability Rate (%)','FontSize',12,'FontName','Arial')
set(gca,'xlim',[0,11])
print(fig,['Figures/','Strain_Viability_Rate'],'-dpng','-r600');  




% Generate figure on Success Rate (Figure 5B)
fig=figure();
set(gcf,'Position',[988   144   353*1.15   216]);
falerate=0.25;
for i=1:length(indexes) %get failure rates 
    fails=yvalues2(:,i)>falerate;
    failures(i)=sum(fails)/length(yvalues2(:,i));
end
bg=barh(1,100*(1-failures(3)));
set(bg,'FaceColor','y')
hold on
bg=barh(2,100*(1-failures(1)));
set(bg,'FaceColor','g')% 
bg=barh(3,100*(1-failures(2)));
set(bg,'FaceColor','b')%
% bg(3).FaceColor='y';
% bg(1).FaceColor='b';
% bg(2).FaceColor='g';
str={['\fontsize{12}C','\fontsize{9}ORI','\fontsize{12}^-'],['\fontsize{12}S','\fontsize{9}LOW'],['\fontsize{12}Q','\fontsize{9}UICK']};
set(gca, 'YtickLabel',str, 'Ytick',1:numel(str),'FontSize',12)
xlabel('Success Rate (%)','FontSize',12,'FontName','Arial')
% set(gca,'xlim',[0.2,3.8])
% set(gca,'ylim',[0 100])
% xtickangle(45)
print(fig,['Figures/','QandSandC Success Rate.png'],'-dpng','-r600');  


end