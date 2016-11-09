%this program estimate the RNAP elongation rate
%based on the transcription initiation delay of three and five prime

%data set information
Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults';

%basic information of data sets
DataInfoFile = fullfile(Folder,filesep,'NewControlDataInfo.xlsx');
[XLSnum5,XLStext5,~] = xlsread(DataInfoFile, 1, '', 'basic');
[XLSnum3,XLStext3,~] = xlsread(DataInfoFile, 2, '', 'basic');

selectedInx5 = strcmp(XLStext5(:,4),'yes');
dataName5 = XLStext5(selectedInx5,1);

selectedInx3 = strcmp(XLStext3(:,4),'yes');
dataName3 = XLStext3(selectedInx3,1);
EffectiveDataSet3 = length(dataName3);

for i0 = 1:length(dataName5);
    FolderNames5{i0} = strrep(dataName5{i0},'\','-');
end

for i0 = 1:length(dataName3);
    FolderNames3{i0} = strrep(dataName3{i0},'\','-');
end

StartTimeFive = cell(3,1);
StartTimeThree = cell(3,1);
for i1 = 1:3
    StartTimeFive{i1} = NaN(41,length(dataName5));
    StartTimeThree{i1} = NaN(41,length(dataName3));
end

temperCountFive = 0;
temperCountThree = 0;


Temperature5 = []; %store all the real temperature for each date set
Temperature3 = []; %store all the real temperature for each date set
for i0 = 1:length(dataName5);
     
        if exist([Folder,filesep,FolderNames5{i0},[filesep,'MeanFitsLoadingRate.mat']])
            load([Folder,filesep,FolderNames5{i0},[filesep,'MeanFitsLoadingRate.mat']]);
            Temperature5 = [Temperature5;XLSnum5(i0)];
            
            temperCountFive = temperCountFive + 1;
            [ROW,COL] = size(FitResults);
            for j0 = 1:ROW
                for k0 = 1:COL
                    if(FitResults(j0,k0).Approved ==1)
                        StartTimeFive{3-COL + k0}(j0,temperCountFive) = FitResults(j0,k0).TimeStart;                       
                    end
                end
            end
        end                   
end


for i0 = 1:EffectiveDataSet3;
     
        if exist([Folder,filesep,FolderNames3{i0},[filesep,'MeanFitsLoadingRate.mat']])
            load([Folder,filesep,FolderNames3{i0},[filesep,'MeanFitsLoadingRate.mat']]);
            Temperature3 = [Temperature3;XLSnum3(i0)];
            
            temperCountThree = temperCountThree + 1;
            [ROW,COL] = size(FitResults);
            for j0 = 1:ROW
                for k0 = 1:COL
                    if(FitResults(j0,k0).Approved ==1)
                        StartTimeThree{3-COL + k0}(j0,temperCountThree) = FitResults(j0,k0).TimeStart;                       
                    end
                end
            end
        end                   
end

UniqueTemperature = unique(Temperature3);
UniqueTempLegend = cell(length(UniqueTemperature),1);
for i0 = 1:length(UniqueTemperature);
    UniqueTempLegend{i0} = num2str(UniqueTemperature(i0));
end

MS2Length = 3.355;

%average over same temperature
MeanStartTimeAP5 = cell(3,1);
StdStartTimeAP5 = cell(3,1);
MeanStartTimeAP3 = cell(3,1);
StdStartTimeAP3 = cell(3,1);
MeanStartTimeFive = nan(length(UniqueTemperature),3);
MeanStartTimeThree = nan(length(UniqueTemperature),3);

SelectedInx = 1:17; %only the anterior part will be counted

for i0 = 1:length(UniqueTemperature);
    SameTempInx = find(Temperature5==UniqueTemperature(i0));
    for j0=1:3;
        MeanStartTimeAP5{j0}(:,i0) = nanmean(StartTimeFive{j0}(SelectedInx,SameTempInx),2);
        StdStartTimeAP5{j0}(:,i0) = nanstd(StartTimeFive{j0}(SelectedInx,SameTempInx),0,2);
        SameTempData = StartTimeFive{j0}(SelectedInx,SameTempInx);
        MeanStartTimeFive(i0,j0) = mean(SameTempData(~isnan(SameTempData)));
    end
end

for i0 = 1:length(UniqueTemperature);
    SameTempInx = find(Temperature3==UniqueTemperature(i0));
    for j0=1:3;
        MeanStartTimeAP3{j0}(:,i0) = nanmean(StartTimeThree{j0}(SelectedInx,SameTempInx),2);
        StdStartTimeAP3{j0}(:,i0) = nanstd(StartTimeThree{j0}(SelectedInx,SameTempInx),0,2);
        SameTempData = StartTimeThree{j0}(SelectedInx,SameTempInx);
        MeanStartTimeThree(i0,j0) = mean(SameTempData(~isnan(SameTempData)));
    end
end


%calculate elongation rates,only considering nc 13 and 14
ElongationRate = cell(2,1);
AverageElongation2 = nan(2,length(UniqueTemperature));
StdAverageElongation2 = nan(2,length(UniqueTemperature));
for i0 = 1:2
    ElongationRate{i0} = nan(size(MeanStartTimeAP5{1},1),size(MeanStartTimeAP5{1},2));
end


%rate along AP 
for k0 = 1:2;
for i0 = 1:size(MeanStartTimeAP5{1},1)
    for j0 = 1:size(MeanStartTimeAP5{1},2)
        if(~isnan(MeanStartTimeAP5{k0+1}(i0,j0)) && ~isnan(MeanStartTimeAP3{k0+1}(i0,j0)))
            ElongationRate{k0}(i0,j0) = MS2Length/(MeanStartTimeAP3{k0+1}(i0,j0)-...
                MeanStartTimeAP5{k0+1}(i0,j0));
        end
    end
end
AverageElongation2(k0,:) = nanmean(ElongationRate{k0},1);
StdAverageElongation2(k0,:) = nanstd(ElongationRate{k0},0,1);
end


TranscriptionDelay = MeanStartTimeThree-MeanStartTimeFive;
AverageElongationRate = MS2Length./TranscriptionDelay(:,2:3);

APAxis = ((0:1:40)*0.025)';
APAxisSelected = (SelectedInx-1)*0.025;
TemperatureLegend = {'18','25.5','30'};

%transcription start time along A-P for three and five prime
%nc 13
figure(1)
subplot(1,2,1)
title('5 Prime','FontSize',20)
hold on
for i0 = 1:size(MeanStartTimeAP5{2},2)
    errorbar(APAxisSelected,MeanStartTimeAP5{2}(:,i0),StdStartTimeAP5{2}(:,i0),'o-','MarkerSize',10,'LineWidth',1)
end
hold off
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'LineWidth',1)
set(gca,'ylim',[2,9],'xlim',[0,0.5])
legend(UniqueTempLegend)

subplot(1,2,2)
title('3 Prime','FontSize',20)
hold on
for i0 = 1:size(MeanStartTimeAP3{2},2)
    errorbar(APAxisSelected,MeanStartTimeAP3{2}(:,i0),StdStartTimeAP3{2}(:,i0),'o-','MarkerSize',10,'LineWidth',1)
end
hold off
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'LineWidth',1)
set(gca,'ylim',[2,9],'xlim',[0,0.5])
legend(UniqueTempLegend)

%nc 14
figure(2)
subplot(1,2,1)
title('5 Prime','FontSize',20)
hold on
for i0 = 1:size(MeanStartTimeAP5{3},2)
    errorbar(APAxisSelected,MeanStartTimeAP5{3}(:,i0),StdStartTimeAP5{3}(:,i0),'o-','MarkerSize',10,'LineWidth',1)
end
hold off
xlabel('AP Axis','FontSize',16,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',20,'LineWidth',1)
set(gca,'ylim',[2,12],'xlim',[0,0.5])
legend(UniqueTempLegend)

subplot(1,2,2)
title('3 Prime','FontSize',20)
hold on
for i0= 1:size(MeanStartTimeAP3{3},2)
    errorbar(APAxisSelected,MeanStartTimeAP3{3}(:,i0),StdStartTimeAP3{3}(:,i0),'o-','MarkerSize',10,'LineWidth',1)
end
hold off
xlabel('AP Axis','FontSize',16,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',20,'LineWidth',1)
set(gca,'ylim',[2,12],'xlim',[0,0.5])
legend(UniqueTempLegend)


%elongation rate along AP axis
%nc 13
figure(3)
title('nc 13','FontSize',24)
plot(APAxisSelected,ElongationRate{1},'o-','MarkerSize',15)
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Elongation rate(kb/min)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',1)
set(gca,'xlim',[0.1,0.6])
legend(TemperatureLegend)

figure(4)
title('nc 14','FontSize',24)
plot(APAxisSelected,ElongationRate{2},'o-','MarkerSize',15)
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Elongation rate(kb/min)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(TemperatureLegend)

%average elongation rates of nc 13 and 14
figure(5)
plot(UniqueTemperature,AverageElongationRate,'o-','LineWidth',2,'MarkerSize',15)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
ylabel('Average elongation rate(kb/min)','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')
legend('nc 13','nc 14')

%average elongation rate, directly average over A-P bins
figure(6)
hold on
errorbar(UniqueTemperature',AverageElongation2(1,:),StdAverageElongation2(1,:),'o-','LineWidth',2,'MarkerSize',15)
errorbar(UniqueTemperature',AverageElongation2(2,:),StdAverageElongation2(2,:),'o-','LineWidth',2,'MarkerSize',15)
hold off
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
ylabel('Average elongation rate(kb/min)','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')
legend('nc 13','nc 14')

