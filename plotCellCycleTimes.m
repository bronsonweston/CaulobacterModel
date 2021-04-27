function [] = plotCellCycleTimes()
% Generate Figure 6A-C of main text
% Weston et al. 2021, Cell Systems

gray=1*[169,169,169]./255;
close all
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)
set(0,'DefaultLegendAutoUpdate','off')

%List of file names containing data.
listOfFiles={'SlowParams_All_SW.mat', 'SlowParams_All_ST.mat', ...
    'CompleteParams_All_SW.mat','CompleteParams_All_ST.mat', ...
    'CtrABindingParams_All_SW.mat', 'CtrABindingParams_All_ST.mat'};
%Allocate memory for importing data
ploty1=zeros(length(listOfFiles)*2,1);
err1=zeros(length(listOfFiles)*2,1);
ploty2=zeros(length(listOfFiles)*2,1);
err2=zeros(length(listOfFiles)*2,1);
ploty3=zeros(length(listOfFiles)*2,1);
err3=zeros(length(listOfFiles)*2,1);

%Import data and data manipulation.
for f=1:length(listOfFiles)
    file=listOfFiles{f};
    load(strcat('ParamCatalog/ChangingBindingSims/',file));
    orig_tchromList=rmoutliers(orig_tchromList); %remove outliers from results of simulation without modification 
    mTcrOri=mean(orig_tchromList); %get mean value
    change_tchromList=rmoutliers(change_tchromList); %remove outliers from results of simulation with CtrAu:Cori binding modification 
    mTcrChange=mean(change_tchromList); %get mean value
    ploty1(2*f-1:2*f)=[mTcrOri, mTcrChange]; %format data into plotable format
    err1(2*(f-1)+1:2*f)=[std(orig_tchromList),std(change_tchromList)]; %get error bar data
end


%% fig 2
for f=1:length(listOfFiles)
    file=listOfFiles{f};
    load(strcat('ParamCatalog/ChangingBindingSims/',file));
    orig_tdivList=rmoutliers(orig_tdivList); %remove outliers from results of simulation without modification 
    mTcrOri=mean(orig_tdivList); %get mean value
    change_tdivList=rmoutliers(change_tdivList);  %remove outliers from results of simulation with CtrAu:Cori binding modification 
    mTcrChange=mean(change_tdivList); %get mean value
    ploty2(2*(f-1)+1:2*f)=[mTcrOri, mTcrChange];  %format data into plotable format
    err2(2*(f-1)+1:2*f)=[std(orig_tchromList),std(change_tchromList)]; %get error bar data
end

%% fig 3
for f=1:length(listOfFiles)
    file=listOfFiles{f};
    load(strcat('ParamCatalog/ChangingBindingSims/',file));
    tchromDifList=rmoutliers(tchromDifList); %remove outliers from results of simulation without modification
    mTcrOri=mean(tchromDifList); %get mean value
    tdivDifList=rmoutliers(tdivDifList); %remove outliers from results of simulation with CtrAu:Cori binding modification 
    mTcrChange=mean(tdivDifList); %get mean value
    ploty3(2*f-1:2*f)=[mTcrOri, mTcrChange]; %format data into plotable format
    err3(2*(f-1)+1:2*f)=[std(tchromDifList),std(tdivDifList)]; %get error bar data
end

%plot graphs
Colors1=[0 0 0; gray; gray]; %define colors of each bar
barChart(ploty1,err1, [],'t^{cr} (min)',Colors1,{'WT','WT +/- CtrA_U:{\itCori}'},[630 215],'CycleTimes1')
barChart(ploty2,err2, [90 180],'t^{div} (min)',Colors1,{'WT','WT +/- CtrA_U:{\itCori}'},[630 215],'CycleTimes2')
barChart(ploty3,err3,[],'\Deltat (min)',Colors1,{'\Deltat^{cr}','\Deltat^{div}'},[630 325],'CycleTimes3')

end

