%this program estimate the RNAP elongation rate

%data set information
Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults';

%basic information of data sets
DataInfoFile = fullfile(Folder,filesep,'NewControlDataInfo.xlsx');
[TemperatureFivePrime,~,~] = xlsread(DataInfoFile, 2, '', 'basic');
[TemperatureThreePrime,~,~] = xlsread(DataInfoFile, 2, '', 'basic');

%load the date
FivePrimeStartTime = load('5primeStartTime.mat');
ThreePrimeStartTime = load('3primeStartTime.mat');

Temperature = unique(TemperatureFivePrime);
MS2Length = 3.355;

%only the anterior part will be counted
SelectedInx = 1:18;

% MeanElongationRate = zero

ElongationRate = cell(2,1);
for i0 = 1:2
    ElongationRate{i0} = nan(size(FivePrimeStartTime.MeanCompareStartTime{1},1),size(FivePrimeStartTime.MeanCompareStartTime{1},2));
end


%rate along AP
for k0 = 1:2;
for i0 = 1:size(FivePrimeStartTime.MeanCompareStartTime{1},1)
    for j0 = 1:size(FivePrimeStartTime.MeanCompareStartTime{1},2)
        if(~isnan(ThreePrimeStartTime.MeanCompareStartTime{k0+1}(i0,j0)) && ~isnan(FivePrimeStartTime.MeanCompareStartTime{k0+1}(i0,j0)))
            ElongationRate{k0}(i0,j0) = MS2Length/(ThreePrimeStartTime.MeanCompareStartTime{k0+1}(i0,j0)-...
                FivePrimeStartTime.MeanCompareStartTime{k0+1}(i0,j0));
        end
    end
end
end


%average elongation rate, average over A-P axis
MeanStartTimeFive = nan(length(Temperature),3);
MeanStartTimeThree = nan(length(Temperature),3);
for j0 = 1:3;
    for i0 = 1:length(Temperature);
        SameTempInxFive = find(TemperatureFivePrime==Temperature(i0));
        SameTempData = FivePrimeStartTime.compareStartTime{j0}(SelectedInx,SameTempInxFive);
        MeanStartTimeFive(i0,j0) = mean(SameTempData(~isnan(SameTempData)));
    end
    
    for i0 = 1:3;
        SameTempInxThree = find(TemperatureThreePrime==Temperature(i0));
        SameTempData = ThreePrimeStartTime.compareStartTime{j0}(SelectedInx,SameTempInxThree);
        MeanStartTimeThree(i0,j0) = mean(SameTempData(~isnan(SameTempData)));
    end
end

TranscriptionDelay = MeanStartTimeThree-MeanStartTimeFive;

APAxis = ((0:1:40)*0.025)';


TemperatureLegend = {'18','25','30'};

%elongation rate along AP axis
figure(1)
plot(APAxis,ElongationRate{1},'o-','MarkerSize',10)
xlabel('AP Axis','FontSize',16,'FontWeight','Bold')
ylabel('Estimated elongation rate(kb/min)','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(TemperatureLegend)

%average elongation rates of nc 13 and 14
figure(2)
plot(Temperature,MS2Length./TranscriptionDelay(:,2:3),'ro','LineWidth',2,'MarkerSize',15)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
ylabel('Average elongation rate(kb/min)','FontSize',16,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])