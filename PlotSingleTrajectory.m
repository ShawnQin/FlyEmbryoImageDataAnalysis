%this program plot single fluorescence trace of the last three data sets
%

clc
clear

%load the data
Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults';
Prefix = '2016-09-15_1-MCP-5-P2P';
load([Folder,filesep,Prefix,filesep,'CompiledParticles.mat']);



figure(1)
hold on
for i0 = 1:length(CompiledParticles)
    plot(ElapsedTime(CompiledParticles(i0).Frame),CompiledParticles(i0).Fluo,'Color', (15-CompiledParticles(i0).nc)/3*[0.8 0.8 0.8])
end
xlabel('Time(min)','FontSize',24,'FontWeight','Bold')
ylabel('Fluorescence','FontSize',24,'FontWeight','Bold')
set(gca,'LineWidth',1.5,'FontSize',24)
hold off