%********************************************************************
%PROGRAM: historicalDataDevTime.m
%DISCRIPTION:
%   this program reproduce some historical data about the temperature
%   dependent development time of various insects
%LAST REVISED: July 19,2016
%********************************************************************
close all
clear
clc

%load the data 
Folder = '/Users/shan/Documents/GoogleDrive/HernanLab/dataAndFigure';
[FileName, DataPath]= uigetfile(Folder,'Choose folder with files to analyze');
DataFile = fullfile(DataPath,FileName);
RawFlyHatch = xlsread(DataFile,1);
RawFlyPupal = xlsread(DataFile,2);
RawEgg = xlsread(DataFile,3);
RawKuniella = xlsread(DataFile,4);

%fitted function 
FitFunc = @(x,K,a,b) K*(1+exp(a+b*x));

%hatch time of fly egg
K1=1/0.070953;
a1=4.4514215;
b1=-0.2071879;
f1 = @(x) K1*(1+exp(a1+b1*x));
xRange1 = (floor(RawFlyHatch(1,1)):0.1:ceil(RawFlyHatch(end,1)))';

%pupal stage of fly
K2=100/36.76;
a2 = 4.435754;
b2 = -0.20036;
f2 = @(x) K2*(1+exp(a2+b2*x));
xRange2 = (floor(RawFlyPupal(1,1)):0.1:ceil(RawFlyPupal(end,1)))';
heatThreshold1 = 29.5;

%incubation time of four different insect
K3 = [100/15.17;100/17.173;100/15.695;100/15.568];
a3 = [4.480468;4.665104;4.308795;4.107601];
b3 = [-0.172223;-0.176068;-0.157424;-0.159571];
xRange3 = (12:0.1:44)';
heatThreshold3 = 34;
FitValue = [];
for i0 = 1:4;
    FitValue = [FitValue,FitFunc(xRange3,K3(i0),a3(i0),b3(i0))];
end

%incubation time of egg Ephestia Kuhniella
K4 = 100/36.35;
a4 = 4.595759;
b4 = -0.206667;
xRange4 = (12:0.1:36)';
heatThreshold4 = 32;
FormulaString4 = '$$y=2.751(1+e^{4.595759-0.206667x})$$';
TextPosi = [25,60];

%test the scaling of hatch time and pupa stage interval hypothesis
%since the two data set are not taken under same temperature, I have to make some
%approximation
SameTemp = [15;20;22.5;25;26;29;30];
newDevTime = nan(length(SameTemp),2);
flag = 1;  %index the Same Temperature
for i0 = 1:length(SameTemp)
    newDevTime(i0,1) = mean(RawFlyHatch(abs(RawFlyHatch(:,1)-SameTemp(i0))<0.4,2));
    newDevTime(i0,2) = newDevTime(i0,1) + mean(RawFlyPupal(abs(RawFlyPupal(:,1)-SameTemp(i0))<0.4,2))*24;
end
RelativeTime = newDevTime(:,1)./newDevTime(:,2);
        

%plot the data

%hatch time of drosophila egg
figure(1)
plot(RawFlyHatch(:,1),RawFlyHatch(:,2),'o','MarkerSize',10)
hold on
plot(xRange1,f1(xRange1),'LineWidth',2)
y1=get(gca,'ylim');
plot([heatThreshold1,heatThreshold1],y1,'k--','LineWidth',1)
% title('nc 12','FontSize',30);
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',1)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',30,'FontWeight','Bold')
ylabel('Hatch Time(hour)','FontSize',30,'FontWeight','Bold')
hold off

%pupal stage time 
figure(2)
plot(RawFlyPupal(:,1),RawFlyPupal(:,2),'o','MarkerSize',10)
hold on
plot(xRange2,f2(xRange2),'LineWidth',2)
y2=get(gca,'ylim');
plot([heatThreshold1,heatThreshold1],y2,'k--','LineWidth',1)
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',1)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',30,'FontWeight','Bold')
ylabel('Pupal stage time (days)','FontSize',30,'FontWeight','Bold')
hold off

%incubation time of four other insects
figure(3)
plot(RawEgg(:,1),RawEgg(:,2:5),'o','MarkerSize',10)
hold on
plot(xRange3,FitValue,'LineWidth',2)
y3=get(gca,'ylim');
plot([heatThreshold3,heatThreshold3],y3,'k--','LineWidth',1)
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',1)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',30,'FontWeight','Bold')
ylabel('Incubation Time(hour)','FontSize',30,'FontWeight','Bold')
hold off


%incubation time of egg Ephestia Kuhniella
figure(4)
plot(RawKuniella(:,1),RawKuniella(:,2),'o','MarkerSize',10)
hold on
plot(xRange4,FitFunc(xRange4,K4,a4,b4),'LineWidth',2)
y4=get(gca,'ylim');
plot([heatThreshold4,heatThreshold4],y4,'k--','LineWidth',1)
text(15,20,FormulaString4,'Interpreter','latex','FontSize',20)
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',1)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',30,'FontWeight','Bold')
ylabel('Incubation Time(hour)','FontSize',30,'FontWeight','Bold')
hold off


%test the scaliing hypothesis
figure(5)
plot(SameTemp,newDevTime,'o','MarkerSize',15,'LineWidth',2)
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',1)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',30,'FontWeight','Bold')
ylabel('Time(hour)','FontSize',30,'FontWeight','Bold')

figure(6)
plot(SameTemp,RelativeTime,'o','MarkerSize',15,'LineWidth',2)
set(gca,'FontSize',24,'FontWeight','Bold','LineWidth',1)
set(gca,'ylim',[0,0.2])
hold on
refline(0,0.17)
hold off
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',30,'FontWeight','Bold')
ylabel('Relative time of hatching','FontSize',30,'FontWeight','Bold')


