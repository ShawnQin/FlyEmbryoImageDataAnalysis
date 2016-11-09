%this program load the fitted RNAP loading rate from single trajectory fitting
%and calculate the average loading rate along AP


close all
clear

Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults';

%load the last three data sets which are taken with higher resolution
Prefix = {'2016-09-16_2-MCP-5-P2P','2016-09-15_1-MCP-5-P2P','2016-09-16_3-MCP-5-P2P'};
MeanData = cell(3,1);
AverageLoading = nan(3,2);
Temperature = [18,25.5,30];

for  k0 = 1:length(Prefix)
load(['/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults',filesep,Prefix{k0},[filesep,'SingleTrajectoryFits.mat']])

ApprovedData = [];
for i0 = 1:length(FitResults)
    if(FitResults(i0).Approved == 1)
        ApprovedData = [ApprovedData;[FitResults(i0).APID,FitResults(i0).TimeStart,FitResults(i0).RateFit,...
            FitResults(i0).ElongationFit,FitResults(i0).nc]];
    end
end

APUnique = unique(ApprovedData(:,1));
MeanData{k0} = zeros(length(APUnique),6);

%only consider nc 14
SlectedIndex = ApprovedData(:,5)==14;
AverageLoading(k0,:) = [mean(ApprovedData(SlectedIndex,3)),std(ApprovedData(SlectedIndex,3))];

for j0 = 1:length(APUnique);
%     APInx = all([ApprovedData(:,1)==APUnique(j0), ApprovedData(:,5)==14],2);
    APInx = ApprovedData(:,1)==APUnique(j0);
    
    
    MeanData{k0}(j0,1:2) = [mean(ApprovedData(APInx,3)),std(ApprovedData(APInx,3))];
    MeanData{k0}(j0,3:4) = [mean(ApprovedData(APInx,4)),std(ApprovedData(APInx,4))];
    MeanData{k0}(j0,5:6) = [mean(ApprovedData(APInx,2)),std(ApprovedData(APInx,2))];
    MeanData{k0}(j0,7) = APUnique(j0)*0.025;
end
end


%plot the results

%loaidng rate along A-P under different temperature
figure(1)
hold on
for i0 = 1:length(MeanData)
    errorbar(MeanData{i0}(:,7),MeanData{i0}(:,1),MeanData{i0}(:,2),'o-','MarkerSize',10,'LineWidth',1.5)   
end 
xlabel('AP Axis','FontSize',24,'FontWeight','Bold')
ylabel('RNAP loading rate(kb/min)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',1)
% set(gca,'xlim',[0.1,0.6])
legend('18','25.5','30')


figure(2)
hold on
for i0 = 1:length(MeanData)
    errorbar(MeanData{i0}(:,7),MeanData{i0}(:,3),MeanData{i0}(:,4),'o-','MarkerSize',10,'LineWidth',1.5)   
end 

figure(3)
errorbar(Temperature,AverageLoading(:,1),AverageLoading(:,2),'ro')
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
ylabel('Average loading rate(kb/min)','FontSize',16,'FontWeight','Bold')