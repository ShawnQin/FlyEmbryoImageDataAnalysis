%**********************************************************
%PROGRAM: compLoadingRateV3.m
%DISCRIPTION:
%       Tis program is a simplified version of "compLoadingRateV2.m", which simply
%       read a list of interested data set, the temperature of each data set has
%       alread specified in "NewControlDataInfo.xlsx'
%       This is used to compare loading rate of new data sets
%LAST REVISED: September 21,2016
%**********************************************************
function compLoadingRateV3(varargin)

close all


%primary folder of the data
Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults';

%basic information of data sets
DataInfoFile = fullfile(Folder,filesep,'NewControlDataInfo.xlsx');
[XLSnum,XLStext,XLSraw] = xlsread(DataInfoFile, 1, '', 'basic');
selectedInx = strcmp(XLStext(:,4),'yes');
dataName = XLStext(selectedInx,1);
EffectiveDataSet = length(dataName);
for i0 = 1:length(dataName);
    FolderNames{i0} = strrep(dataName{i0},'\','-');
end
% Temperature = unique(XLSnum);  %unique temperature  


compareLoadingRate = cell(3,1); %for different nuclei cycle

%initialization
for i1 = 1:3
    compareLoadingRate{i1} = NaN(41,EffectiveDataSet);
end

temperCount = 0;



Temperature = []; %store all the real temperature for each date set
for i0 = 1:EffectiveDataSet
     
        if exist([Folder,filesep,FolderNames{i0},[filesep,'MeanFitsLoadingRate.mat']])
            load([Folder,filesep,FolderNames{i0},[filesep,'MeanFitsLoadingRate.mat']]);
            Temperature = [Temperature;XLSnum(i0)];
            
            %initialize
            %only consider those bins approved
            temperCount = temperCount + 1;
            [ROW,COL] = size(FitResults);
            for j0 = 1:ROW
                for k0 = 1:COL
                    if(FitResults(j0,k0).Approved ==1)
                        compareLoadingRate{3-COL + k0}(j0,temperCount) = FitResults(j0,k0).RateFit;                       
                    end
                end
            end
        end
               
        
end

%mean value for data sets taken at same temperature

%legend
[UniqueTemperature,IndexTemp,OrderTemp] =  unique(Temperature);
UniqueTempLegend = cell(length(UniqueTemperature),1);
for i0 = 1:length(UniqueTemperature);
    UniqueTempLegend{i0} = num2str(UniqueTemperature(i0));
end

% NumberRepeat = [UniqueTemperature(2:end),]
MeanCompareLoadingRate = cell(3,1);
StdCompareLoadingRate = cell(3,1);
SummaryLoadingRate = cell(3,1);

if ~isempty(varargin)
    if(strcmp(varargin{1},'mid'))
        SelectedInx = 15:41;
    elseif(strcmp(varargin{1},'anterior'))
        SelectedInx = 1:16; 
    elseif(strcmp(varargin{1},'max'))
        SelectedInx = 5:16;
    end
else
    SelectedInx = 1:41;  %all the AP 
end

for i0 = 1:length(UniqueTemperature);
    SameTempInx = find(Temperature==UniqueTemperature(i0));
    for j0=1:3;
%         nanmean(compareLoadingRate{j0}(:,SameTempInx),2)
        MeanCompareLoadingRate{j0}(:,i0) = nanmean(compareLoadingRate{j0}(SelectedInx,SameTempInx),2);
        StdCompareLoadingRate{j0}(:,i0) = nanstd(compareLoadingRate{j0}(SelectedInx,SameTempInx),0,2);
        SummaryLoadingRate{j0}(i0,1) = nanmean(reshape(compareLoadingRate{j0}(SelectedInx,SameTempInx),1,[]));
        SummaryLoadingRate{j0}(i0,2) = nanstd(reshape(compareLoadingRate{j0}(SelectedInx,SameTempInx),1,[]));
        
    end
end

%plot the results
%loading rate under different temperature
figure(1)
hold on
APAxis = (SelectedInx*0.025)';
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
% APAxis = (0:1:40)'*0.025;
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
% APAxis = (0:1:40)'*0.025;
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






%mean loading rate over AP
figure(4)
%nc 12
% subplot(2,2,1)
errorbar(UniqueTemperature',SummaryLoadingRate{1}(:,1),SummaryLoadingRate{1}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 12','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Relative loading rate of RNAP','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
hold on

%nc 13
% subplot(2,2,2)
errorbar(UniqueTemperature',SummaryLoadingRate{2}(:,1),SummaryLoadingRate{2}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 13','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Relative loading rate of RNAP','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%nc 14
% subplot(2,2,3)
errorbar(UniqueTemperature',SummaryLoadingRate{3}(:,1),SummaryLoadingRate{3}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 14','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Relative loading rate of RNAP','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
legend('nc 12','nc 13','nc 14')
hold off
