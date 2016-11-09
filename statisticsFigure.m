%****************************************************
%PRGRAM NAME: statisticsFigure.m
%DESCRIPTION: the data is a structure called CompiliedParticles.mat
%             modified from the 2015 KITP summer course
%LAST REVISED: Aug 10,2016
%
%*********************************************************

function statisticsFigure(varargin)



%load the data
Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults';
FolderTemp = uigetdir(Folder,'Choose folder with files to analyze');
if exist([FolderTemp,filesep,'CompiledParticles.mat'])
   load([FolderTemp,filesep,'CompiledParticles.mat']);
end


%plot the mean fluoresence as a function of time with error bar
figure(1)
% errorbar((1:length(MeanVectorAll))/3,MeanVectorAll,SDVectorAll)
hold on
% plot((1:length(MeanVectorAll))/3,MeanVectorAll,'o-','MarkerSize',8,'LineWidth',2)

ParticleNumber = length(CompiledParticles);
for j = 1:ParticleNumber;
    if CompiledParticles(j).nc==13
    plot(ElapsedTime(CompiledParticles(j).Frame),CompiledParticles(j).Fluo,'color',[0.8 0.8 0.8])
    end
end

plot(ElapsedTime,MeanVectorAll,'-','MarkerSize',8,'LineWidth',3)


title('temperature 20.26','FontSize',24)
xlabel('time(min)','FontSize',24,'FontWeight','Bold')
ylabel('fluorescence(a.u)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
hold off



%plot single particles
% totFrame = length(MeanVectorAll);
% ParticleNumber = length(CompiledParticles);
% allSingleFluo = cell(ParticleNumber,1);
% frames = cell(ParticleNumber,1);
% newSingleFluo = cell(ParticleNumber,1);
% for i = 1:ParticleNumber;
%     allSingleFluo{i} = CompiledParticles(i).Fluo;
%     frames{i} = CompiledParticles(i).Frame;
%     newSingleFluo{i} = nan(1,totFrame);
%     newSingleFluo{i}(frames{i})=allSingleFluo{i};
% end
% % 
% %single particle profile
% figure(2)
% hold on
% for j = 1:ParticleNumber;
%     plot(ElapsedTime(frames{j}),CompiledParticles(j).Fluo,'color',[0.8 0.8 0.8])
% end
% xlabel('time(min)','FontSize',24,'FontWeight','Bold')
% ylabel('fluorecence(a.u)','FontSize',24,'FontWeight','Bold')
% set(gca,'FontSize',16,'FontWeight','Bold')
% hold off


Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data';

DataInfoFile = fullfile(Folder,filesep,'DynamicsResults',filesep,'ResultsInfo.xlsx');
[XLSnum,XLStext,XLSraw] = xlsread(DataInfoFile, 1, '', 'basic');
AllDate = datenum(datetime(XLSnum(:,1),'ConvertFrom','excel1900'));

%load the data for analysis
FolderTemp = uigetdir(Folder,'Choose folder with files to analyze');
AllFolder = dir(FolderTemp);
[uniqueC,~,idx] = unique(XLStext(:,4));
counts = accumarray(idx(:),1,[],@sum);
EffectiveDataSet = max(counts);  %this is the number of useful data sets

for j0 = 1:length(AllFolder) 
    FolderNames{j0} = AllFolder(j0).name;
end 

%initialization
compareMaximumIntensity = cell(3,1);
compareStartTime = cell(3,1);  %transcription start time after mitois
compareTotalmRNA = cell(3,1); %for different nuclei cycle

for i0 = 1:3;
    compareMaximumIntensity{i0} = NaN(41,EffectiveDataSet);
    compareStartTime{i0} = NaN(41,EffectiveDataSet);
    compareTotalmRNA{i0} = NaN(41,EffectiveDataSet);
end

Temperature = [];
temperCount = 0;
for i0 = 1:length(FolderNames)
    
    if(isdir([FolderTemp,filesep,FolderNames{i0}]))
        %comparing loading rate
        DataFolder = [FolderTemp,filesep,AllFolder(i0).name];
        if(exist([DataFolder,[filesep,'CompiledParticles.mat']])...
                && exist([DataFolder,[filesep,'APDivision.mat']]))           
            [MeanMaxIntensityAP,MeanStartTimeAP,SelectedTrace,meanmRNA]=...
                maxIntenStartTime(DataFolder);
            
            temperCount = temperCount + 1;
            for j0 = 1:3;
                compareMaximumIntensity{j0}(:,temperCount) = MeanMaxIntensityAP(1+j0,:)';
                compareStartTime{j0}(:,temperCount) = MeanStartTimeAP(1+j0,:)';
                compareTotalmRNA{j0}(:,temperCount) = meanmRNA(1+j0,:)';
            end
            %find out the temperature of this data set
            HypenPosition = find(FolderNames{i0}=='-');
            dateString = FolderNames{i0}(1:HypenPosition(3)-1);
            dateOfThisSet = datenum(dateString,'yyyy-mm-dd');
            INX = find(AllDate==dateOfThisSet);
            Temperature = [Temperature;XLSnum(INX,2)];         
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
MeanMaximumIntensity = cell(3,1);
StdMaximumIntensity = cell(3,1);

MeanStartTime = cell(3,1);
StdStartTime = cell(3,1);


MeanCompareTotalmRNA = cell(3,1);
StdCompareTotalmRNA = cell(3,1);

SummaryMaximumIntensity = cell(3,1);
SummarymRNA = cell(3,1);
SummaryStartTime = cell(3,1);


if ~isempty(varargin)
    if(strcmp(varargin{1},'mid'))
        SelectedInx = 15:41;
    elseif(strcmp(varargin{1},'anterior'))
        SelectedInx = 1:16; 
    end
else
    SelectedInx = 1:41;  %all the AP 
end

for i0 = 1:length(UniqueTemperature);
    SameTempInx = find(Temperature==UniqueTemperature(i0));
    for j0=1:3;
%         nanmean(compareLoadingRate{j0}(:,SameTempInx),2)
        MeanMaximumIntensity{j0}(:,i0) = nanmean(compareMaximumIntensity{j0}(SelectedInx,SameTempInx),2);
        StdMaximumIntensity{j0}(:,i0) = nanstd(compareMaximumIntensity{j0}(SelectedInx,SameTempInx),0,2);
        SummaryMaximumIntensity{j0}(i0,1) = nanmean(reshape(compareMaximumIntensity{j0}(SelectedInx,SameTempInx),1,[]));
        SummaryMaximumIntensity{j0}(i0,2) = nanstd(reshape(compareMaximumIntensity{j0}(SelectedInx,SameTempInx),1,[]));
        
        MeanStartTime{j0}(:,i0) = nanmean(compareStartTime{j0}(SelectedInx,SameTempInx),2);
        StdStartTime{j0}(:,i0) = nanstd(compareStartTime{j0}(SelectedInx,SameTempInx),0,2);
        SummaryStartTime{j0}(i0,1) = nanmean(reshape(compareStartTime{j0}(SelectedInx,SameTempInx),1,[]));
        SummaryStartTime{j0}(i0,2) = nanstd(reshape(compareStartTime{j0}(SelectedInx,SameTempInx),1,[]));
       
        MeanCompareTotalmRNA{j0}(:,i0) = nanmean(compareTotalmRNA{j0}(SelectedInx,SameTempInx),2);
        StdCompareTotalmRNA{j0}(:,i0) = nanstd(compareTotalmRNA{j0}(SelectedInx,SameTempInx),0,2);
        SummarymRNA{j0}(i0,1) = nanmean(reshape(compareTotalmRNA{j0}(SelectedInx,SameTempInx),1,[]));
        SummarymRNA{j0}(i0,2) = nanstd(reshape(compareTotalmRNA{j0}(SelectedInx,SameTempInx),1,[]));
        
    end
end



%comparing the maximum intensity of fluorescence
%nc 12
APAxis = (SelectedInx*0.025)';
figure(1)
hold on
for i0 = 1:size(MeanMaximumIntensity{1},2)
    errorbar(APAxis,MeanMaximumIntensity{1}(:,i0),StdMaximumIntensity{1}(:,i0),'o-','MarkerSize',...
        10,'LineWidth',1)
end
hold off
title('nc 12','FontSize',30);
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('maximum fluorscence intensity(a.u)','FontSize',24,'FontWeight','Bold')
legend(UniqueTempLegend)

%nc 13
figure(2)
hold on
for i0 = 1:size(MeanMaximumIntensity{2},2)
    errorbar(APAxis,MeanMaximumIntensity{2}(:,i0),StdMaximumIntensity{2}(:,i0),'o-','MarkerSize',...
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
figure(3)
hold on
for i0 = 1:size(MeanMaximumIntensity{3},2)
    errorbar(APAxis,MeanMaximumIntensity{3}(:,i0),StdMaximumIntensity{3}(:,i0),'o-','MarkerSize',...
        10,'LineWidth',1)
end
hold off
title('nc 14','FontSize',30);
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',24,'FontWeight','Bold')
legend(UniqueTempLegend)





%comparing transcription start time 
%nc 12
figure(4)
hold on
for i0 = 1:size(MeanStartTime{1},2)
    errorbar(APAxis,MeanStartTime{1}(:,i0),StdStartTime{1}(:,i0),'o-','MarkerSize',...
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
figure(5)
hold on
for i0 = 1:size(MeanStartTime{2},2)
    errorbar(APAxis,MeanStartTime{2}(:,i0),StdStartTime{2}(:,i0),'o-','MarkerSize',...
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
figure(6)
hold on
for i0 = 1:size(MeanStartTime{3},2)
    errorbar(APAxis,MeanStartTime{3}(:,i0),StdStartTime{3}(:,i0),'o-','MarkerSize',...
        10,'LineWidth',1)
end
hold off
title('nc 14','FontSize',30);
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',24,'FontWeight','Bold')
legend(UniqueTempLegend)



%plot the fitted total mRNA produced during cycle 12 and 13
figure(7)
%nc 12
subplot(2,2,1)
hold on
for i0 = 1:size(MeanCompareTotalmRNA{1},2)
    errorbar(APAxis,MeanCompareTotalmRNA{1}(:,i0),StdCompareTotalmRNA{1}(:,i0))
end
hold off
title('nc 12','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Total mRNA (a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)

%nc 13
subplot(2,2,2)
hold on
for i0 = 1:size(MeanCompareTotalmRNA{2},2)
    errorbar(APAxis,MeanCompareTotalmRNA{2}(:,i0),StdCompareTotalmRNA{2}(:,i0))
end
hold off
legend()
title('nc 13','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Total mRNA (a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)


%nc 14
subplot(2,2,3)
hold on
for i0 = 1:size(MeanCompareTotalmRNA{3},2)
    errorbar(APAxis,MeanCompareTotalmRNA{3}(:,i0),StdCompareTotalmRNA{3}(:,i0))
end
hold off
legend()
title('nc 13','FontSize',20)
xlabel('AP Axis','FontSize',20,'FontWeight','Bold')
ylabel('Total mRNA (a.u)','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')
set(gca,'xlim',[0.1,0.6])
legend(UniqueTempLegend)



%summary of maximum intensity, start time and mRNA accumulation
figure(8)
%average maximum intensity
subplot(2,2,1)
hold on
for i0 = 1:3
    errorbar(UniqueTemperature,SummaryMaximumIntensity{i0}(:,1),SummaryMaximumIntensity{i0}(:,2),'o-')
end
hold off
title('average maximum intensity','FontSize',20)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',16,'FontWeight','Bold')
ylabel('fluorescence intensity (a.u)','FontSize',24,'FontWeight','Bold')


subplot(2,2,2)
%average transcription start time
hold on
for i0 = 1:3
    errorbar(UniqueTemperature,SummaryStartTime{i0}(:,1),SummaryStartTime{i0}(:,2),'o-')
end
hold off
title('average transcription start time ','FontSize',20)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',16,'FontWeight','Bold')
ylabel('Start time(min)','FontSize',24,'FontWeight','Bold')

subplot(2,2,3)
%average total mRNA
hold on
for i0 = 1:3
    errorbar(UniqueTemperature,SummarymRNA{i0}(:,1),SummarymRNA{i0}(:,2),'o-')
end
hold off
title('average mRNA accumulation','FontSize',20)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',16,'FontWeight','Bold')
ylabel('mRNA accumulated(a.u)','FontSize',24,'FontWeight','Bold')

end


