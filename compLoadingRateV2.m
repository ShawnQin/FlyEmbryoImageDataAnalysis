%**********************************************************
%PROGRAM: compLoadingRateV2.m
%DISCRIPTION:
%       Tis program is a simplified version of "compLoadingRate.m", which simply
%       read a list of interested data set, the temperature of each data set has
%       alread specified in "NewControlDataInfo.xlsx'
%       This is used to compare loading rate of new data sets
%LAST REVISED: September 5th,2016
%**********************************************************
function compLoadingRateV2(varargin)

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
compareStartTime = cell(3,1);  %transcription start time after mitois
compareTotalmRNA = cell(3,1);  %for different nuclei cycle
compareElongationRate = cell(3,1); %estimated elongation rate for different nc 
compareOnTime = cell(2,1);    %estimated ON time during cycle 12 and 13
compareOffRate = cell(2,1);   %no off rate of nuclei cycle 14
compareFitTotalmRNA = cell(2,1);  %the total mRNA estimated are based on the fitted parameters

% information about the construct
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
for i0 = 1:EffectiveDataSet
     
        if exist([Folder,filesep,FolderNames{i0},[filesep,'MeanFitsLoadingElongationRate.mat']])
            load([Folder,filesep,FolderNames{i0},[filesep,'MeanFitsLoadingElongationRate.mat']]);
            Temperature = [Temperature;XLSnum(i0)];
            
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
        if exist([Folder,filesep,FolderNames{i0},[filesep,'MeanmRNA.mat']])
           load([Folder,filesep,FolderNames{i0},[filesep,'MeanmRNA.mat']]);
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
MeanCompareStartTime = cell(3,1);
StdCompareStartTime = cell(3,1);
MeanCompareTotalmRNA = cell(3,1);
StdCompareTotalmRNA = cell(3,1);
MeanCompareElongationRate = cell(3,1);
StdCompareElongationRate = cell(3,1);

MeanCompareOffRate = cell(2,1); %Off rate of nc 12 and 13
StdCompareOffRate = cell(2,1);
MeanCompareFitTotalmRNA = cell(2,1);  %fited total mRNA produced during nc 12 and 13
StdCompareFitTotalmRNA = cell(2,1);
MeanCompareOnTime = cell(2,1);   %fitted ON time during nc 12 and 13
StdCompareOnTime = cell(2,1); 

SummaryLoadingRate = cell(3,1);
SummarymRNA = cell(3,1);
SummaryElongationRate = cell(3,1);
SummaryStartTime = cell(3,1);
SummaryOnTime = cell(3,1);

if ~isempty(varargin)
    if(strcmp(varargin{1},'mid'))
        SelectedInx = 15:41;
    elseif(strcmp(varargin{1},'anterior'))
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
        
        MeanCompareStartTime{j0}(:,i0) = nanmean(compareStartTime{j0}(SelectedInx,SameTempInx),2);
        StdCompareStartTime{j0}(:,i0) = nanstd(compareStartTime{j0}(SelectedInx,SameTempInx),0,2);
        SummaryStartTime{j0}(i0,1) = nanmean(reshape(compareStartTime{j0}(SelectedInx,SameTempInx),1,[]));
        SummaryStartTime{j0}(i0,2) = nanstd(reshape(compareStartTime{j0}(SelectedInx,SameTempInx),1,[]));
       
        MeanCompareTotalmRNA{j0}(:,i0) = nanmean(compareTotalmRNA{j0}(SelectedInx,SameTempInx),2);
        StdCompareTotalmRNA{j0}(:,i0) = nanstd(compareTotalmRNA{j0}(SelectedInx,SameTempInx),0,2);
        SummarymRNA{j0}(i0,1) = nanmean(reshape(compareTotalmRNA{j0}(SelectedInx,SameTempInx),1,[]));
        SummarymRNA{j0}(i0,2) = nanstd(reshape(compareTotalmRNA{j0}(SelectedInx,SameTempInx),1,[]));
        
        MeanCompareElongationRate{j0}(:,i0) = nanmean(compareElongationRate{j0}(SelectedInx,SameTempInx),2);
        StdCompareElongationRate{j0}(:,i0) = nanstd(compareElongationRate{j0}(SelectedInx,SameTempInx),0,2);
        SummaryElongationRate{j0}(i0,1) = nanmean(reshape(compareElongationRate{j0}(SelectedInx,SameTempInx),1,[]));
        SummaryElongationRate{j0}(i0,2) = nanstd(reshape(compareElongationRate{j0}(SelectedInx,SameTempInx),1,[]));
%     MeanCompareLoadingRate{j0}(:,i0) = compareLoadingRate{j0}(:,SameTempInx);
%     MeanCompareTotalmRNA{j0}(:,i0)  = compareTotalmRNA{j0}(:,compareTotalmRNA);
    end
    
    for k0 = 1:2
        MeanCompareFitTotalmRNA{k0}(:,i0) = nanmean(compareFitTotalmRNA{k0}(SelectedInx,SameTempInx),2);
        StdCompareFitTotalmRNA{k0}(:,i0) = nanstd(compareFitTotalmRNA{k0}(SelectedInx,SameTempInx),0,2);
        
        MeanCompareOnTime{k0}(:,i0) = nanmean(compareOnTime{k0}(SelectedInx,SameTempInx),2);
        StdCompareOnTime{k0}(:,i0) = nanstd(compareOnTime{k0}(SelectedInx,SameTempInx),0,2);
        SummaryOnTime{k0}(i0,1) = nanmean(reshape(compareOnTime{k0}(SelectedInx,SameTempInx),1,[]));
        SummaryOnTime{k0}(i0,2) = nanstd(reshape(compareOnTime{k0}(SelectedInx,SameTempInx),1,[]));
        
        MeanCompareOffRate{k0}(:,i0) = nanmean(compareOffRate{k0}(SelectedInx,SameTempInx),2);
        StdCompareOffRate{k0}(:,i0) = nanstd(compareOffRate{k0}(SelectedInx,SameTempInx),0,2);
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
legend(UniqueTempLegend)

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
legend(UniqueTempLegend)

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
legend(UniqueTempLegend)






%comparing transcription start time 
%nc 12
figure(8)
hold on
for i0 = 1:size(MeanCompareStartTime{1},2)
    errorbar(APAxis,MeanCompareStartTime{1}(:,i0),StdCompareStartTime{1}(:,i0),'o-','MarkerSize',...
        10,'LineWidth',1)
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
        10,'LineWidth',1)
end
hold off
title('nc 13','FontSize',30);
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',24,'FontWeight','Bold')
legend(UniqueTempLegend)

%nc 14
figure(10)
hold on
for i0 = 1:size(MeanCompareStartTime{1},2)
    errorbar(APAxis,MeanCompareStartTime{3}(:,i0),StdCompareStartTime{3}(:,i0),'o-','MarkerSize',...
        10,'LineWidth',1)
end
hold off
title('nc 14','FontSize',30);
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',24,'FontWeight','Bold')
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
ylabel('Duration of transcription(min)','FontSize',20,'FontWeight','Bold')
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
ylabel('Duration of transcription(min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)




%plot elongation rate 
figure(13)
%nc 12
subplot(2,2,1)
hold on
for i0 = 1:size(MeanCompareElongationRate{1},2)
    errorbar(APAxis,MeanCompareElongationRate{1}(:,i0),StdCompareElongationRate{1}(:,i0))
end
hold off
title('nc 12','FontSize',20)
xlabel('AP Axis','FontSize',16,'FontWeight','Bold')
ylabel('Estimated elongation rate','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)

%nc 13
subplot(2,2,2)
hold on
for i0 = 1:size(MeanCompareElongationRate{2},2)
    errorbar(APAxis,MeanCompareElongationRate{2}(:,i0),StdCompareElongationRate{2}(:,i0),'o-','MarkerSize',10)
end
hold off
title('nc 13','FontSize',20)
xlabel('AP Axis','FontSize',16,'FontWeight','Bold')
ylabel('Estimated elongation rate(kb/min)','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)

%nc 14
subplot(2,2,3)
hold on
for i0 = 1:size(MeanCompareElongationRate{3},2)
    errorbar(APAxis,MeanCompareElongationRate{3}(:,i0),StdCompareElongationRate{3}(:,i0),'o-','MarkerSize',10)
end
hold off
title('nc 14','FontSize',20)
xlabel('AP Axis','FontSize',16,'FontWeight','Bold')
ylabel('Estimated elongation rate(kb/min)','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)




%plot Off rate
figure(14)
%nc 12
subplot(1,2,1)
hold on
for i0 = 1:size(MeanCompareOffRate{1},2)
    errorbar(APAxis,abs(MeanCompareOffRate{1}(:,i0)),StdCompareOffRate{1}(:,i0))
end
hold off
title('nc 12','FontSize',20)
xlabel('AP Axis','FontSize',16,'FontWeight','Bold')
ylabel('Off rate(a.u)','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)

%nc 13
subplot(1,2,2)
hold on
for i0 = 1:size(MeanCompareOffRate{2},2)
    errorbar(APAxis,abs(MeanCompareOffRate{2}(:,i0)),StdCompareOffRate{2}(:,i0))
end
hold off
title('nc 13','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Off rate','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)




%mean loading rate over AP
figure(17)
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



%plot the overall elongation rate, sum over AP 
figure(18)
%nc 12
% subplot(2,2,1)
errorbar(UniqueTemperature',SummaryElongationRate{1}(:,1),SummaryElongationRate{1}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 12','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Elongation rate of RNAP (kb/min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
hold on

%nc 13
% subplot(2,2,2)
errorbar(UniqueTemperature',SummaryElongationRate{2}(:,1),SummaryElongationRate{2}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 13','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Elongation rate of RNAP (kb/min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%nc 14
% subplot(2,2,3)
errorbar(UniqueTemperature',SummaryElongationRate{3}(:,1),SummaryElongationRate{3}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 14','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Elongation rate of RNAP (kb/min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
legend('nc 12','nc 13','nc 14')
hold off




%plot overall start time
figure(19)
%nc 12
subplot(2,2,1)
errorbar(UniqueTemperature',SummaryStartTime{1}(:,1),SummaryStartTime{1}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 12','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Transcription initiation time(min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%nc 13
subplot(2,2,2)
errorbar(UniqueTemperature',SummaryStartTime{2}(:,1),SummaryStartTime{2}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 13','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Transcription initiation time(min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%nc 14
subplot(2,2,3)
errorbar(UniqueTemperature',SummaryStartTime{3}(:,1),SummaryStartTime{3}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 14','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Transcription initiation time(min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')




%plot overall On time
figure(20)
%nc 12
% subplot(1,2,1)
errorbar(UniqueTemperature',SummaryOnTime{1}(:,1),SummaryOnTime{1}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 12','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Duration of On state(min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
hold on

%nc 13
% subplot(1,2,2)
errorbar(UniqueTemperature',SummaryOnTime{2}(:,1),SummaryOnTime{2}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 13','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Duration of On state(min)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
hold off





%plot the mean overall mRNA accumulated during each NC
figure(21)
%nc 12
subplot(2,2,1)
errorbar(UniqueTemperature',SummarymRNA{1}(:,1),SummarymRNA{1}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 12','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Mean mRNA acculumation(a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%nc 13
subplot(2,2,2)
errorbar(UniqueTemperature',SummarymRNA{2}(:,1),SummarymRNA{2}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 13','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Mean mRNA acculumation(a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

%nc 14
subplot(2,2,3)
errorbar(UniqueTemperature',SummarymRNA{3}(:,1),SummarymRNA{3}(:,2),'o-','MarkerSize',15,'LineWidth',2)
title('nc 14','FontSize',20)
xlabel('Temperature','FontSize',20,'FontWeight','Bold')
ylabel('Mean mRNA acculumation(a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')

