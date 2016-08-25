%PROGRAM: CompareDivisionTimeMaxIntensity.m
%DISCRIPTION: This program compare the division time along A-P axis as well as the
%     maximum fluorescence intensity, which represent the spacing of RNAP on DNA
%     load date from 'CompiledParticles.mat' and 'APDivision.mat'
%     The exact duration of each nc is then calculated from the 'ElapsedTime.mat'
%LAST REVISED: Aug 7,2016

close
clear
clc

%load the data for analysis
%primary folder of the data
Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data';
FolderTemp = uigetdir(Folder,'Choose folder with files to analyze');
AllFolder = dir(FolderTemp);


%temperature information of each of each data set
DataInfoFile = fullfile(Folder,filesep,'DynamicsResults',filesep,'ResultsInfo.xlsx');
[XLSnum,XLStext,XLSraw] = xlsread(DataInfoFile, 1, '', 'basic');
AllDate = datenum(datetime(XLSnum(:,1),'ConvertFrom','excel1900'));
% [FileName, DataPath]= uigetfile(Folder,'Choose folder with files to analyze');
% DataInfoFile = fullfile(DataPath,FileName);
% [XLSnum,XLStex,XLSraw] = xlsread(DataInfoFile, 1, '', 'basic');
% AllDate = datenum(datetime(XLSnum(:,1),'ConvertFrom','excel1900'));

ncLength = [];
Temperature = [];
temperCount2 = 0;
for i0 = 1:length(AllFolder)
    if(isdir([FolderTemp,filesep,AllFolder(i0).name]))
        if (exist([FolderTemp,filesep,AllFolder(i0).name,[filesep,'CompiledParticles.mat']])...
                && exist([FolderTemp,filesep,AllFolder(i0).name,[filesep,'APDivision.mat']]))
           load([FolderTemp,filesep,AllFolder(i0).name,[filesep,'CompiledParticles.mat']]);
           load([FolderTemp,filesep,AllFolder(i0).name,[filesep,'APDivision.mat']]);
           
           %find the temperature of this date set
           HypenPosition = find(AllFolder(i0).name=='-');
           dateString = AllFolder(i0).name(1:HypenPosition(3)-1);
           dateOfThisSet = datenum(dateString,'yyyy-mm-dd');
           INX = find(AllDate==dateOfThisSet);
           Temperature = [Temperature;XLSnum(INX,2)];
           temperCount2 = temperCount2 +1;
           
%            newAPDivision = nan(size(APDivision,1),size(APDivision,2));
           newAPDivision = APDivision;
           newAPDivision(newAPDivision~=0) = ElapsedTime(APDivision(APDivision~=0));
           newAPDivision(newAPDivision==0) = nan;
%            newAPDivision(APDivision~=0) = APDivision(APDivision~=0);
           for i3=1:3;
               ncLength(temperCount2,i3)= nanmean(newAPDivision(11+i3,:)-newAPDivision(10+i3,:));
           end
        end
    end
end

[UniqueTemperature,IndexTemp,OrderTemp] =  unique(Temperature);
UniqueEstiTemp = 10.6556 + 0.5651*UniqueTemperature;
UniqueTempLegend = [];
for i0 = 1:length(UniqueTemperature);
    UniqueTempLegend = [UniqueTempLegend;num2str(UniqueEstiTemp(i0))];
end
MeanNCLength = nan(length(UniqueTemperature),3);
StdNCLength = nan(length(UniqueTemperature),3);
for i0 = 1:length(UniqueTemperature);
    SameTempInx = find(Temperature==UniqueTemperature(i0));
    MeanNCLength(i0,:) = nanmean(ncLength(SameTempInx,:),1);
    StdNCLength(i0,:) = nanstd(ncLength(SameTempInx,:),1);
end

%proportion of each nc length
TotLength = sum(ncLength,2);
StdLen = [];
MeanLen = nanmean(ncLength,1);
for i0 = 1:3;
    RatioTime(:,i0) = ncLength(~isnan(TotLength),i0)./TotLength(~isnan(TotLength));
    StdLen = [StdLen,nanstd(ncLength(:,i0))];
end
Proportion = cumsum(RatioTime,2);
TemperatureThreeNC = Temperature(~isnan(TotLength));

%random generize the number with the same statistics of ncLength
NumRand = 12;
RandLength = zeros(NumRand,3);
for i0 = 1:3
    RandLength(:,i0) = normrnd(MeanLen(i0),StdLen(i0),NumRand,1);
end

RandRaio = cumsum(RandLength,2);
for i0 = 1:3
    RandGenerateProportion(:,i0) = RandRaio(:,i0)./sum(RandLength,2);
end


%plot the nc length as a function of temperature
%cycle 13
figure(1)
errorbar(UniqueEstiTemp,MeanNCLength(:,2),StdNCLength(:,2),'o','MarkerSize',10,'LineWidth',2)
title('nc 13','FontSize',30);
set(gca,'FontSize',24,'FontWeight','Bold')
% set(gca,'xlim',[0.1,0.6])
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
ylabel('Nuclei Cycle Length(min)','FontSize',24,'FontWeight','Bold')
% legend(UniqueTempLegend)

%cycle 14
figure(2)
errorbar(UniqueEstiTemp,MeanNCLength(:,3),StdNCLength(:,3),'o','MarkerSize',10,'LineWidth',2)
title('nc 14','FontSize',30);
set(gca,'FontSize',24,'FontWeight','Bold')
% set(gca,'xlim',[0.1,0.6])
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
ylabel('Nuclei Cycle Length(min)','FontSize',24,'FontWeight','Bold')
% legend(UniqueTempLegend)



%scatter length nc 12 and 13
figure(3)
subplot(1,2,1)
hold on
plot(ncLength(:,1),ncLength(:,2),'o','MarkerSize',12,'LineWidth',2)
plot(RandLength(:,1),RandLength(:,2),'o','MarkerSize',12,'LineWidth',2)
hold off
xlabel('nc lenght of 12 (min)','FontSize',24,'FontWeight','Bold')
ylabel('nc length of 13 (min)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
legend('experiment','random')

subplot(1,2,2)
hold on
plot(ncLength(:,2),ncLength(:,3),'o','MarkerSize',12,'LineWidth',2)
plot(RandLength(:,2),RandLength(:,3),'o','MarkerSize',12,'LineWidth',2)
hold off
xlabel('nc lenght of 13 (min)','FontSize',24,'FontWeight','Bold')
ylabel('nc length of 14 (min)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
legend('experiment','random')




%test the proportion of each nc length
figure(4)
hold on
for i0 = 1:length(TemperatureThreeNC);
    plot(Proportion(i0,:),ones(1,3)*TemperatureThreeNC(i0),'o','MarkerSize',12,'LineWidth',2)
end
hold off
xlabel('Proportion of last three cycles lengths','FontSize',24,'FontWeight','Bold')
ylabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')


%generated random data
figure(6)
hold on 
for i0  = 1:NumRand
    plot(RandGenerateProportion(i0,:),ones(1,3)*i0,'o','MarkerSize',12,'LineWidth',2)
end
hold off
xlabel('Proportion of last three cycles lengths','FontSize',24,'FontWeight','Bold')
ylabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')







