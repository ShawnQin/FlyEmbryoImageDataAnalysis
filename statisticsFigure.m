%****************************************************
%PRGRAM NAME: statisticsFigure.m
%DESCRIPTION: the data is a structure called CompiliedParticles.mat
%             modified from the 2015 KITP summer course
%LAST REVISED: July 10,2016
%
%*********************************************************
clf
clear
clc

%load the data
Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults';
for i = 1:2;  %compare two data set taken at same temperature
FolderTemp = uigetdir(Folder,'Choose folder with files to analyze');
if exist([FolderTemp,filesep,'CompiledParticles.mat'])
   load([FolderTemp,filesep,'CompiledParticles.mat']);
end


%plot the mean fluoresence as a function of time with error bar
figure(1)
% errorbar((1:length(MeanVectorAll))/3,MeanVectorAll,SDVectorAll)
hold on
plot((1:length(MeanVectorAll))/3,MeanVectorAll,'o-','MarkerSize',8,'LineWidth',2)
end
hold off
title('temperature 20.26','FontSize',24)
xlabel('time(min)','FontSize',24,'FontWeight','Bold')
ylabel('fluorescence(a.u)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',20,'FontWeight','Bold')



%plot single particles
totFrame = length(MeanVectorAll);
ParticleNumber = 101;
allSingleFluo = cell(ParticleNumber,1);
frames = cell(ParticleNumber,1);
newSingleFluo = cell(ParticleNumber,1);
for i = 1:ParticleNumber;
    allSingleFluo{i} = CompiledParticles(i).Fluo;
    frames{i} = CompiledParticles(i).Frame;
    newSingleFluo{i} = nan(1,totFrame);
    newSingleFluo{i}(frames{i})=allSingleFluo{i};
end

%single particle profile
figure(2)
hold on
particleLable = []; %store the particles that last more than 100 frames
for j = 1:ParticleNumber;
    if(length(allSingleFluo{j})>=100)
        plot((1:totFrame)/6,newSingleFluo{j})
        particleLable = [particleLable;j];
    end
end
xlabel('time(min)','FontSize',24,'FontWeight','Bold')
ylabel('fluorecence(a.u)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')

hold off



