%PROGRAM: CompareDivisionTimeMaxIntensity.m
%DISCRIPTION: This program compare the division time along A-P axis as well as the
%     maximum fluorescence intensity, which represent the spacing of RNAP on DNA
%     load date from 'CompiledParticles.mat' and 'APDivision.mat'
%     The exact duration of each nc is then calculated from the 'ElapsedTime.mat'
%     Revision 1: add arguments at the beginning, can either assign some data sets
%     or import the whole folder
%LAST REVISED: Sept 5,2016

function CompareDivisionTime(varargin)

%load the data for analysis
%primary folder of the data
Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data';
ncLength = [];

if ~isempty(varargin)
    if strcmp(varargin,'specify')
       FolderName = {'2016-09-03_1-MCP-5-P2P','2016-09-03_2-MCP-5-P2P','2016-09-06-MCP-5-P2P','2016-09-09_2-MCP-5-P2P',...
        '2016-09-10_2-MCP-5-P2P','2016-09-10_3-MCP-5-P2P','2016-09-13-MCP-5-P2P','2016-09-14_1-MCP-3-P2P',...
        '2016-09-14_2-MCP-3-P2P','2016-09-13-MCP-5-P2P'};
       Temperature = [25.5,25.5,30,25.5,18,30,30];
       
    elseif strcmp(varargin,'load')
        DataInfoFile = fullfile(Folder,filesep,'DynamicsResults',filesep,'NewControlDataInfo.xlsx');
        [XLSnum,XLStext,XLSraw] = xlsread(DataInfoFile,3, '', 'basic');
        dataName = XLStext(2:end,1);
        Temperature = XLSnum;
        EffectiveDataSet = length(dataName);
        for i0 = 1:length(dataName);
            FolderName{i0} = strrep(dataName{i0},'\','-');
        end
    end
    
    UniqueEstiTemp = unique(Temperature);
    for i0 = 1:length(FolderName)
        if(exist([Folder,filesep,'DynamicsResults',filesep,FolderName{i0},[filesep,'CompiledParticles.mat']])...
                && exist([Folder,filesep,'DynamicsResults',filesep,FolderName{i0},[filesep,'APDivision.mat']]))
              load([Folder,filesep,'DynamicsResults',filesep,FolderName{i0},[filesep,'CompiledParticles.mat']]);
              load([Folder,filesep,'DynamicsResults',filesep,FolderName{i0},[filesep,'APDivision.mat']]);
              newAPDivision = APDivision;
              newAPDivision(newAPDivision~=0) = ElapsedTime(APDivision(APDivision~=0));
              newAPDivision(newAPDivision==0) = nan;
            for i3=1:3;
               ncLength(i0,i3)= nanmean(newAPDivision(11+i3,:)-newAPDivision(10+i3,:));
            end
        end
    end
    
        

else
    %temperature information of each of each data set
    DataInfoFile = fullfile(Folder,filesep,'DynamicsResults',filesep,'ResultsInfo.xlsx');
    [XLSnum,XLStext,XLSraw] = xlsread(DataInfoFile, 1, '', 'basic');
    AllDate = datenum(datetime(XLSnum(:,1),'ConvertFrom','excel1900'));

    FolderTemp = uigetdir(Folder,'Choose folder with files to analyze');
    AllFolder = dir(FolderTemp);
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
end


UniqueTempLegend = cell(length(UniqueEstiTemp),1);
for i0 = 1:length(UniqueEstiTemp);
    UniqueTempLegend{i0} = num2str(UniqueEstiTemp(i0));
end
MeanNCLength = nan(length(UniqueEstiTemp),3);
StdNCLength = nan(length(UniqueEstiTemp),3);
for i0 = 1:length(UniqueEstiTemp);
    SameTempInx = find(Temperature==UniqueEstiTemp(i0));
    MeanNCLength(i0,:) = nanmean(ncLength(SameTempInx,:),1);
    StdNCLength(i0,:) = nanstd(ncLength(SameTempInx,:),1,1);
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
% NumRand = 12;
% RandLength = zeros(NumRand,3);
% for i0 = 1:3
%     RandLength(:,i0) = normrnd(MeanLen(i0),StdLen(i0),NumRand,1);
% end
% 
% RandRaio = cumsum(RandLength,2);
% for i0 = 1:3
%     RandGenerateProportion(:,i0) = RandRaio(:,i0)./sum(RandLength,2);
% end

%theoretical extrapolation of nc length
ScalingFun = @(T1,T2,alpha) exp(alpha./T1)/exp(alpha./T2);
alpha = 37.31;
Temperature = (15:2:33);
RefTemp = 25;
RefNCLength = [9.7,12.6,16]; %data from Foe 1989 paper


%based on earlier study
%fitted function 
a1=4.4514215;
b1=-0.2071879;
FitFunc = @(T1,T2,a,b) (1+exp(a+b*T1))./(1+exp(a+b*T2));

allNCLength = zeros(length(Temperature),3);
NCLength2 = zeros(length(Temperature),3);
for i0 = 1:length(Temperature);
    allNCLength(i0,:) = ScalingFun(Temperature(i0),RefTemp,alpha).*RefNCLength;
    NCLength2(i0,:) = FitFunc(Temperature(i0),RefTemp,a1,b1).*RefNCLength;
end



%plot the nc length as a function of temperature
%cycle 13
figure(1)
hold on
errorbar(UniqueEstiTemp,MeanNCLength(:,2),StdNCLength(:,2),'o','MarkerSize',10,'LineWidth',2)
title('nc 12','FontSize',30);
plot(Temperature',allNCLength(:,2),'k--','LineWidth',2)
plot(Temperature',NCLength2(:,2),'k-.','LineWidth',2)
set(gca,'FontSize',24,'FontWeight','Bold')
% set(gca,'xlim',[0.1,0.6])
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
ylabel('Nucleus Cycle Length(min)','FontSize',24,'FontWeight','Bold')
legend('experiment','prediction 1','prediction 2')
hold off
% legend(UniqueTempLegend)

%cycle 14
figure(2)
hold on
errorbar(UniqueEstiTemp,MeanNCLength(:,3),StdNCLength(:,3),'o','MarkerSize',10,'LineWidth',2)
plot(Temperature',allNCLength(:,3),'k--','LineWidth',2)
plot(Temperature',NCLength2(:,3),'k-.','LineWidth',2)
title('nc 13','FontSize',30);
set(gca,'FontSize',24,'FontWeight','Bold')
% set(gca,'xlim',[0.1,0.6])
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
ylabel('Nucleus Cycle Length(min)','FontSize',24,'FontWeight','Bold')
legend('experiment','prediction 1','prediction 2')
hold off
% legend(UniqueTempLegend)



% %scatter length nc 12 and 13
% figure(3)
% subplot(1,2,1)
% hold on
% plot(ncLength(:,1),ncLength(:,2),'o','MarkerSize',12,'LineWidth',2)
% plot(RandLength(:,1),RandLength(:,2),'o','MarkerSize',12,'LineWidth',2)
% hold off
% xlabel('nc lenght of 12 (min)','FontSize',24,'FontWeight','Bold')
% ylabel('nc length of 13 (min)','FontSize',24,'FontWeight','Bold')
% set(gca,'FontSize',20,'FontWeight','Bold')
% legend('experiment','random')
% 
% subplot(1,2,2)
% hold on
% plot(ncLength(:,2),ncLength(:,3),'o','MarkerSize',12,'LineWidth',2)
% plot(RandLength(:,2),RandLength(:,3),'o','MarkerSize',12,'LineWidth',2)
% hold off
% xlabel('nc lenght of 13 (min)','FontSize',24,'FontWeight','Bold')
% ylabel('nc length of 14 (min)','FontSize',24,'FontWeight','Bold')
% set(gca,'FontSize',20,'FontWeight','Bold')
% legend('experiment','random')




% %test the proportion of each nc length
% figure(4)
% hold on
% for i0 = 1:length(TemperatureThreeNC);
%     plot(Proportion(i0,:),ones(1,3)*TemperatureThreeNC(i0),'o','MarkerSize',12,'LineWidth',2)
% end
% hold off
% xlabel('Proportion of last three cycles lengths','FontSize',24,'FontWeight','Bold')
% ylabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
% set(gca,'FontSize',24,'FontWeight','Bold')
% 
% 
% %generated random data
% figure(6)
% hold on 
% for i0  = 1:NumRand
%     plot(RandGenerateProportion(i0,:),ones(1,3)*i0,'o','MarkerSize',12,'LineWidth',2)
% end
% hold off
% xlabel('Proportion of last three cycles lengths','FontSize',24,'FontWeight','Bold')
% ylabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
% set(gca,'FontSize',24,'FontWeight','Bold')

end






