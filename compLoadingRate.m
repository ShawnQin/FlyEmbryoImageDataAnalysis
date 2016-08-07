%**********************************************************
%PROGRAM: compLoadingRate.m
%DISCRIPTION:
%       This program read the data of loading rates under different temperature
%       and comparing temperature dependent.
%       We specify the A-P position, which is 0.025 of the whole embryos in length
%       The method used to estimate cummulated mRNA has revised, the area under the
%       curve is not the correct way to do
%LAST REVISED: August 5,2016
%**********************************************************
function compLoadingRate(varargin)

close all


%primary folder of the data
Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data';

DataInfoFile = fullfile(Folder,filesep,'DynamicsResults',filesep,'ResultsInfo.xlsx');
[XLSnum,XLStext,XLSraw] = xlsread(DataInfoFile, 1, '', 'basic');
AllDate = datenum(datetime(XLSnum(:,1),'ConvertFrom','excel1900'));

if ~isempty(varargin)
    %selection type of data set, by default  it will load all the subfolder in chosen
    %folder, use 'mid' to open only defined folders
    if(strcmp(varargin{1},'mid'))

%   this only consider the middle part of the embryo, we are trying to see the
%   expression boundary of Hb. Need to added the folder name each time you have new
%   data sets. This is a compromised way to do the job.
    FolderNames = {'2016-07-23-MCP-5-P2P','2016-07-24-MCP-5-P2P','2016-07-25-MCP-5-P2P','2016-07-26-MCP-5-P2P',...
        '2016-07-26-MCP-5-P2P','2016-07-27-MCP-5-P2P','2016-07-28-MCP-5-P2P'};
    
    elseif(strcmp(varargin{1},'anterior'))
    FolderNames = {'2016-07-02-MCP-5-P2P','2016-07-05-MCP-5-P2P','2016-07-07-MCP-5-P2P','2016-07-11-MCP-5-P2P',...
        '2016-07-12-MCP-5-P2P','2016-07-14-MCP-5-P2P','2016-07-17-MCP-5-P2P','2016-07-18-MCP-5-P2P'};    
    end
    FolderTemp = [Folder,filesep,'DynamicsResults'];
    EffectiveDataSet = length(FolderNames);
else
    
%load the data for analysis
FolderTemp = uigetdir(Folder,'Choose folder with files to analyze');
AllFolder = dir(FolderTemp);
[uniqueC,~,idx] = unique(XLStext(:,4));
counts = accumarray(idx(:),1,[],@sum);
EffectiveDataSet = max(counts);  %this is the number of useful data sets

for j0 = 1:length(AllFolder) 
    FolderNames{j0} = AllFolder(j0).name;
end 
end

%information file of each data set
%load the data 
% Folder = '/Users/shan/Documents/GoogleDrive/HernanLab/dataAndFigure';
% [FileName, DataPath]= uigetfile(Folder,'Choose folder with files to analyze');



%loop through different temperatures
compareLoadingRate = cell(3,1); %for different nuclei cycle
compareStartTime = cell(3,1);  %transcription start time after mitois
compareTotalmRNA = cell(3,1); %for different nuclei cycle
compareElongationRate = cell(3,1); %estimated elongation rate for different nc 
compareOnTime = cell(2,1); %estimated ON time during cycle 12 and 13
compareOffRate = cell(2,1);   %no off rate of nuclei cycle 14
compareFitTotalmRNA = cell(2,1);  %the total mRNA estimated are based on the fitted parameters


TemperatureSet = [23,17,27,38,30,30,17,23,38];  %different number of temperature
RealTempEstim = 10.6556 + 0.5651*TemperatureSet;
diffTemperature = length(TemperatureSet);   %number of different temperatures tried

MS2Length = 1336;
LacZLength = 3960;
alpha = MS2Length/LacZLength;
%initialization
for i1 = 1:3
    compareLoadingRate{i1} = NaN(41,EffectiveDataSet);
    compareStartTime{i1} = NaN(41,EffectiveDataSet);
    compareTotalmRNA{i1} = NaN(41,EffectiveDataSet);
    compareElongationRate{i1} = NaN(41,EffectiveDataSet);
%     ncMeanLength{i1} = nan(41,diffTemperature);
end

for i0  = 1:2;
    compareOffRate{i0} = NaN(41,EffectiveDataSet);
    compareFitTotalmRNA{i0} = NaN(41,EffectiveDataSet);
    compareOnTime{i0} = NaN(41,EffectiveDataSet);
end

ncMeanLength = [];  %length of nc

temperCount = 0;
temperCount2 = 0;
temperCount3 = 0;


Temperature = []; %store all the real temperature for each date set
for i0 = 1:length(FolderNames)
    
    if(isdir([FolderTemp,filesep,FolderNames{i0}]))
        %comparing loading rate
        if exist([FolderTemp,filesep,FolderNames{i0},[filesep,'MeanFitsLoadingElongationRate.mat']])
            load([FolderTemp,filesep,FolderNames{i0},[filesep,'MeanFitsLoadingElongationRate.mat']]);
            %find out the temperature of this data set
            HypenPosition = find(FolderNames{i0}=='-');
            dateString = FolderNames{i0}(1:HypenPosition(3)-1);
            dateOfThisSet = datenum(dateString,'yyyy-mm-dd');
            INX = find(AllDate==dateOfThisSet);
            Temperature = [Temperature;XLSnum(INX,2)];
            
            %initialize
            
            %only consider those bins approved
            temperCount = temperCount + 1;
            [ROW,COL] = size(FitResults);
            for j0 = 1:ROW
                for k0 = 1:COL
                    if(FitResults(j0,k0).Approved ==1)
                        compareLoadingRate{3-COL + k0}(j0,temperCount) = FitResults(j0,k0).RateFit;
                        compareStartTime{3-COL + k0}(j0,temperCount) = FitResults(j0,k0).TimeStart;
                        compareElongationRate{3-COL + k0}(j0,temperCount) = FitResults(j0,k0).ElongationFit;
                        if(k0<3) %cycle 14 doesn't have off rate
                            compareOffRate{3-COL + k0}(j0,temperCount) = FitResults(j0,k0).RateOffFit;
                            
                            compareOnTime{3-COL + k0}(j0,temperCount) = FitResults(j0,k0).TimeEnd-FitResults(j0,k0).TimeStart;
                            
                            DelayTime = (MS2Length+LacZLength)/1000/FitResults(j0,k0).ElongationFit;
                            PleteauIntensity = DelayTime*FitResults(j0,k0).RateFit;
                            compareFitTotalmRNA{3-COL + k0}(j0,temperCount) = (1+alpha)/(2+alpha)*(compareOnTime{3-COL + k0}(j0,temperCount))/DelayTime*2*PleteauIntensity;
                                        
                        end
                    end
                end
            end
        end
               
        
        %comparing the total mRNA along AP axis
        if exist([FolderTemp,filesep,FolderNames{i0},[filesep,'MeanmRNA.mat']])
           load([FolderTemp,filesep,FolderNames{i0},[filesep,'MeanmRNA.mat']]);
           temperCount3 = temperCount3 +1;
            [ROW,COL] = size(FitResults);
            for j3 = 1:ROW
                for k3 = 1:COL
                    if(FitResults(j3,k3).Approved ==1)
                        try
                            compareTotalmRNA{3-COL + k3}(j3,temperCount3) = FitResults(j3,k3).TotalmRNA;
%                       compareOffRate{3-COL + k0}(j0,temperCount) = FitResults(j0,k0).RateOffFit;
                        catch
%                             compareTotalmRNA{3-COL + k3}(j3,temperCount3) = nan;
                            error('The dimension of the rhs and lhs do not match!')
                        end
                    end
                end
            end
        end
               
           
    end 
    
    
end


%mean value for data sets taken at same temperature

%legend
[UniqueTemperature,IndexTemp,OrderTemp] =  unique(Temperature);
UniqueEstiTemp = 10.6556 + 0.5651*UniqueTemperature;
UniqueTempLegend = [];
for i0 = 1:length(UniqueTemperature);
    UniqueTempLegend = [UniqueTempLegend;num2str(UniqueEstiTemp(i0))];
end
% NumberRepeat = [UniqueTemperature(2:end),]
MeanCompareLoadingRate = cell(3,1);
StdCompareLoadingRate = cell(3,1);
MeanCompareStartTime = cell(3,1);
StdCompareStartTime = cell(3,1);
MeanCompareTotalmRNA = cell(3,1);
StdCompareTotalmRNA = cell(3,1);
MeanCompareElongationRate = cell(3,1);
StdCompareElongationRate = cell(3,1);
MeanCompareFitTotalmRNA = cell(2,1);  %fited total mRNA produced during nc 12 and 13
StdCompareFitTotalmRNA = cell(2,1);
MeanCompareOnTime = cell(2,1);   %fitted ON time during nc 12 and 13
StdCompareOnTime = cell(2,1); 

for i0 = 1:length(UniqueTemperature);
    SameTempInx = find(Temperature==UniqueTemperature(i0));
    for j0=1:3;
%         nanmean(compareLoadingRate{j0}(:,SameTempInx),2)
        MeanCompareLoadingRate{j0}(:,i0) = nanmean(compareLoadingRate{j0}(:,SameTempInx),2);
        StdCompareLoadingRate{j0}(:,i0) = nanstd(compareLoadingRate{j0}(:,SameTempInx),0,2);
        
        MeanCompareStartTime{j0}(:,i0) = nanmean(compareStartTime{j0}(:,SameTempInx),2);
        StdCompareStartTime{j0}(:,i0) = nanstd(compareStartTime{j0}(:,SameTempInx),0,2);

        MeanCompareTotalmRNA{j0}(:,i0) = nanmean(compareTotalmRNA{j0}(:,SameTempInx),2);
        StdCompareTotalmRNA{j0}(:,i0) = nanstd(compareTotalmRNA{j0}(:,SameTempInx),0,2);
        
        MeanCompareElongationRate{j0}(:,i0) = nanmean(compareElongationRate{j0}(:,SameTempInx),2);
        StdCompareElongationRate{j0}(:,i0) = nanstd(compareElongationRate{j0}(:,SameTempInx),0,2);
%     MeanCompareLoadingRate{j0}(:,i0) = compareLoadingRate{j0}(:,SameTempInx);
%     MeanCompareTotalmRNA{j0}(:,i0)  = compareTotalmRNA{j0}(:,compareTotalmRNA);
    end
    
    for k0 = 1:2
        MeanCompareFitTotalmRNA{k0}(:,i0) = nanmean(compareFitTotalmRNA{k0}(:,SameTempInx),2);
        StdCompareFitTotalmRNA{k0}(:,i0) = nanstd(compareFitTotalmRNA{k0}(:,SameTempInx),0,2);
        
        MeanCompareOnTime{k0}(:,i0) = nanmean(compareOnTime{k0}(:,SameTempInx),2);
        StdCompareOnTime{k0}(:,i0) = nanstd(compareOnTime{k0}(:,SameTempInx),0,2);
    end

end

%plot the results
%loading rate under different temperature
figure(1)
hold on
%cycle 12
% subplot(2,2,1)
APAxis = (0:1:40)'*0.025;
% plot(APAxis,compareLoadingRate{1},'o-','MarkerSize',10,'LineWidth',2)
for i0 = 1:size(MeanCompareLoadingRate{1},2)
    errorbar(APAxis,MeanCompareLoadingRate{1}(:,i0),StdCompareLoadingRate{1}(:,i0),'o-','MarkerSize',10,'LineWidth',2)
end
hold off
title('nc 12','FontSize',30);
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Relative loading rate of RNAP','FontSize',24,'FontWeight','Bold')
% legend(TempLegend)
legend(UniqueTempLegend)

%cycle 13
% subplot(2,2,2)
figure(2)
hold on
APAxis = (0:1:40)'*0.025;
% plot(APAxis,compareLoadingRate{2},'o-','MarkerSize',10,'LineWidth',2)
for i0 = 1:size(MeanCompareLoadingRate{2},2)
    errorbar(APAxis,MeanCompareLoadingRate{2}(:,i0),StdCompareLoadingRate{2}(:,i0),'o-','MarkerSize',10,'LineWidth',2)
end
hold off
title('nc 13','FontSize',30);
set(gca,'xlim',[0.1,0.6])
set(gca,'FontSize',20,'FontWeight','Bold')
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Relative loading rate of RNAP','FontSize',24,'FontWeight','Bold')
% legend(TempLegend)
legend(UniqueTempLegend)



%cycle 14
% subplot(2,2,3)
figure(3)
hold on
APAxis = (0:1:40)'*0.025;
% plot(APAxis,compareLoadingRate{3},'o-','MarkerSize',10,'LineWidth',2)
for i0 = 1:size(MeanCompareLoadingRate{3},2)
    errorbar(APAxis,MeanCompareLoadingRate{3}(:,i0),StdCompareLoadingRate{3}(:,i0),'o-',...
        'MarkerSize',10,'LineWidth',2)
end
hold off
title('nc 14','FontSize',30);
set(gca,'xlim',[0.1,0.6])
set(gca,'FontSize',20,'FontWeight','Bold')
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Relative loading rate of RNAP','FontSize',24,'FontWeight','Bold')
% legend(TempLegend)
legend(UniqueTempLegend)


%plot the temperature dependent rate at a certain AP position
APAxisSelect = input('Select a position(0~1):');
% P_inx = round(POSITION/0.025);

% APAxisSelect = 0.35;
whichBin = round(APAxisSelect/0.025);
[TempOrder,Inx] = sort(RealTempEstim);

figure(4)
subplot(1,2,1)
errorbar(UniqueEstiTemp',MeanCompareLoadingRate{2}(whichBin,:),StdCompareLoadingRate{2}(whichBin,:),...
    'o','MarkerSize',8,'LineWidth',1)
title('nc 13','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Relative loading rate of RNAP','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

subplot(1,2,2)
% plot(TempOrder,compareLoadingRate{3}(whichBin,Inx),'o','MarkerSize',10,'LineWidth',2)
errorbar(UniqueEstiTemp',MeanCompareLoadingRate{3}(whichBin,:),StdCompareLoadingRate{3}(whichBin,:),...
    'o','MarkerSize',8,'LineWidth',1)
title('nc 14','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Relative loading rate of RNAP','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%plot the total mRNA along AP Axis during different nc
figure(6)
subplot(2,2,1)
%nc 12
hold on
for i0 = 1:size(MeanCompareTotalmRNA{1},2)
    errorbar(APAxis,MeanCompareTotalmRNA{1}(:,i0),StdCompareTotalmRNA{1}(:,i0))
end
hold off
title('nc 12','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Area under the curve','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])

%nc 13
subplot(2,2,2)
% plot(APAxis,compareTotalmRNA{2},'o-','MarkerSize',8,'LineWidth',2)
hold on
for i0 = 1:size(MeanCompareTotalmRNA{2},2)
    errorbar(APAxis,MeanCompareTotalmRNA{2}(:,i0),StdCompareTotalmRNA{2}(:,i0))
end
hold off
title('nc 13','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Area under the curve','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])

%nc 14
subplot(2,2,3)
% plot(APAxis,compareTotalmRNA{3},'o-','MarkerSize',8,'LineWidth',2)
hold on
for i0 = 1:size(MeanCompareTotalmRNA{3},2)
    errorbar(APAxis,MeanCompareTotalmRNA{3}(:,i0),StdCompareTotalmRNA{3}(:,i0))
end
hold off
title('nc 14','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Area under the curve','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])

%compare at a certain position
figure(7)
POSITION = input('Select a position(0~1):');
P_inx = round(POSITION/0.025);

%nc 12
subplot(2,2,1)
% plot(TempOrder,compareTotalmRNA{1}(P_inx,Inx),'o','MarkerSize',10,'LineWidth',2)
errorbar(UniqueEstiTemp',MeanCompareTotalmRNA{1}(P_inx,:),StdCompareTotalmRNA{1}(P_inx,:),'o',...
    'MarkerSize',10,'LineWidth',2)
title('nc 12','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Area under the curve','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%nc 13
subplot(2,2,2)
% plot(TempOrder,compareTotalmRNA{2}(P_inx,Inx),'o','MarkerSize',10,'LineWidth',2)
errorbar(UniqueEstiTemp',MeanCompareTotalmRNA{2}(P_inx,:),StdCompareTotalmRNA{2}(P_inx,:),'o',...
    'MarkerSize',10,'LineWidth',2)
title('nc 13','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Area under the curve','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%nc 14
subplot(2,2,3)
% plot(TempOrder,compareTotalmRNA{3}(P_inx,Inx),'o','MarkerSize',10,'LineWidth',2)
errorbar(UniqueEstiTemp',MeanCompareTotalmRNA{3}(P_inx,:),StdCompareTotalmRNA{3}(P_inx,:),'o',...
    'MarkerSize',10,'LineWidth',2)
title('nc 14','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Area under the curve','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%nc 12
figure(8)
hold on
for i0 = 1:size(MeanCompareStartTime{1},2)
    errorbar(APAxis,MeanCompareStartTime{1}(:,i0),StdCompareStartTime{1}(:,i0),'o-','MarkerSize',...
        10,'LineWidth',2)
end
hold off
title('nc 12','FontSize',30);
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',24,'FontWeight','Bold')
legend(UniqueTempLegend)

%nc 13
figure(9)
hold on
for i0 = 1:size(MeanCompareStartTime{2},2)
    errorbar(APAxis,MeanCompareStartTime{2}(:,i0),StdCompareStartTime{2}(:,i0),'o-','MarkerSize',...
        10,'LineWidth',2)
end
hold off
title('nc 13','FontSize',30);
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',24,'FontWeight','Bold')
legend(UniqueTempLegend)


figure(10)
hold on
for i0 = 1:size(MeanCompareStartTime{1},2)
    errorbar(APAxis,MeanCompareStartTime{3}(:,i0),StdCompareStartTime{3}(:,i0),'o-','MarkerSize',...
        10,'LineWidth',2)
end
hold off
title('nc 14','FontSize',30);
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',24,'FontWeight','Bold')
legend(UniqueTempLegend)

%plot the fitted total mRNA produced during cycle 12 and 13
figure(11)
%nc 12
subplot(1,2,1)
hold on
for i0 = 1:size(MeanCompareFitTotalmRNA{1},2)
    errorbar(APAxis,MeanCompareFitTotalmRNA{1}(:,i0),StdCompareFitTotalmRNA{1}(:,i0))
end
hold off
title('nc 12','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Estimated total mRNA (a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)

%nc 13
subplot(1,2,2)
hold on
for i0 = 1:size(MeanCompareFitTotalmRNA{2},2)
    errorbar(APAxis,MeanCompareFitTotalmRNA{2}(:,i0),StdCompareFitTotalmRNA{1}(:,i0))
end
hold off
legend()
title('nc 13','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Estimated total mRNA (a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)


%plot the ON time during each nc
figure(12)
%nc 12
subplot(1,2,1)
hold on
for i0 = 1:size(MeanCompareOnTime{1},2)
    errorbar(APAxis,MeanCompareOnTime{1}(:,i0),StdCompareOnTime{1}(:,i0))
end
hold off
title('nc 13','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Estimated total mRNA (a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)

%nc 13
subplot(1,2,2)
hold on
for i0 = 1:size(MeanCompareOnTime{2},2)
    errorbar(APAxis,MeanCompareOnTime{2}(:,i0),StdCompareOnTime{1}(:,i0))
end
hold off
title('nc 13','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Estimated total mRNA (a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)




%save all the data extract for future plot
FoderSaveData = '/Users/shan/Documents/MATLAB/LivemRNAFISH/FinalDataAndFigure';
save([FoderSaveData,filesep,'AllCompRates.mat'],'compareLoadingRate',...
    'compareLoadingRate','compareOffRate','compareOffRate','MeanCompareLoadingRate',...
    'StdCompareLoadingRate','MeanCompareStartTime','StdCompareStartTime','MeanCompareTotalmRNA',...
    'StdCompareTotalmRNA')
