%this program extrapolate nc length from Kuntz 2014 paper
close all
clear

ScalingFun = @(T1,T2,alpha) exp(alpha./T1)/exp(alpha./T2);
alpha = 37.31;
BaseTemp = 273;
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



figure(1)
plot(Temperature,allNCLength,'o-','MarkerSize',12,'LineWidth',2)
hold on
plot(Temperature,NCLength2,'^--','MarkerSize',12,'LineWidth',2)
xlabel(['Temperature ','(',sprintf('%c', char(176)),'C)'],'FontSize',24,'FontWeight','Bold')
ylabel('nc length (min)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')
legend('nc 12','nc 13','nc 14')