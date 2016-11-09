%this program try to analysis single particle trajecoties
%I tried two method: one is aligning all the trajectory to a "average" trace and 
%then using previous fitting method, the other one is directly fit the a two-state
%model to a single trajectory, and then average along A-P bins

%Last revised on 09/30/2016


function SingleTarjectoryFit(varargin)

%modified from 'FitMeanAPLoadingRate2.m', manualy choose the range to do linear fitting
%use changeable elongation rate
%The final results are named 'MeanFitsLoadingElongationRate.mat'

%V2, separate the rate of turning on and rate of turning off. I found that
%the off slope is higher than on suggesting mitotic repression.
%This code will grab results from V1 and convert them to V2 in the first
%round. The mat files will not be backwards compatible.
%I also modified the errors to report the 68% CI.

%This function performs fits to the mean fluorescence as a function of time
%of a particular dataset.
%OUTPUT: MeanFitsAsymmetric.mat
%It gives you n columns each representing a nuclear cycle and m rows each
%representing a bin number

%Fitting:
%a,z: On time
%s,x: Off time
%d,c: Rate, fine
%D,C: Rate, coarse
%F,V: Off Rate, coarse
%f,v: Off Rate, fine
%q,w: approve or disapprove a fit
%e: save
%G,B: elongation rate, coarse
%g,b: elongation rate, fine

%Moving around:
%, .: Move in AP
%n,m: Move in nc
%k,l: Change fit range from the right
%h,j: Change fit range from the left


%Parameters:
% % MinParticles=2;     %Minimum number of particles in an AP bin
MinTimePoints=5;    %Minimum number of time points where we'll have at least
                    %the minimum number of particles.
ElongationRate=1.54;    %In kb/minutes, this is just a reference elongation rate
GeneLength=5.296;       %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb.
Delay=GeneLength/ElongationRate;    %Minutes for PolII to fall off after reaching
                                    %the first MS2 site.

                                    
close all


%Get the default folders
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
%     DetermineLocalFolders;
% 
% if ~isempty(varargin)
%     Prefix=varargin{1};
%                
% else
%     FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
%     Dashes=strfind(FolderTemp,filesep);
%     Prefix=FolderTemp((Dashes(end)+1):end);
% end

%only used when debugging
DropboxFolder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults';
Prefix = '2016-09-16_2-MCP-5-P2P';
load(['/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults',filesep,Prefix,[filesep,'APDivision.mat']])
load(['/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults',filesep,Prefix,[filesep,'CompiledParticles.mat']])


% %Get the relevant folders now:
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
%     DetermineLocalFolders(Prefix);
%         


%Load the complied particles and the division information
%to make this program works for windows, mac and linux systems
% load([DropboxFolder,filesep,Prefix,[filesep,'CompiledParticles.mat']])
% %number of approved particles
NumParticles = length(CompiledParticles);
% 
% if exist([DropboxFolder,filesep,Prefix,[filesep,'APDivision.mat']])
%     load([DropboxFolder,filesep,Prefix,[filesep,'APDivision.mat']])
% else
%     error('Could not load APDivision.mat. Make sure to have done the manual check of division.')
% end




%Initial parameters for fits. We will estimate the maximum rate based on
%the elongation time and the maximum average fluorescence of the data set.
MaxRate=max(max(MeanVectorAP))/Delay;


Rate014=MaxRate;     %Rate per minute
Elongation014 = ElongationRate;
TimeStart014=5;
TimeEnd014=1000;                                    
                           


 
%Rough frame window to consider in the fits

%Some data sets won't have nc12
if nc12>0
    FrameWindow12=[nc12:nc13];
else
    FrameWindow12=[];
end

%some data won't have the meitosis 12
if nc13>0
    FrameWindow13=[nc13:nc14];
else
    FrameWindow13 = (1:nc14);
end
FrameWindow14=[nc14:length(ElapsedTime)];      

         

%Set the first guess for the parameters for each particles and also
%dissaproved the ones that did not have enough data points. Fit results has
%is a structure with the fits corresponding to each particle in AP position and nc13
%or nc14
if exist([DropboxFolder,filesep,Prefix,filesep,'SingleTrajectoryFits.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'SingleTrajectoryFits.mat']);
    
    if isempty(FitResults)
        FitResults(NumParticles,1).Rate0=[];
        FitResults(NumParticles,1).Elongation0=[];
    elseif(~isfield(FitResults,'Elongation0'))
        FitResults(NumParticles,1).Elongation0=[];
    end
else
    FitResults(NumParticles,1).Rate0=[];
    FitResults(NumParticles,1).Elongation0=[];
end




%Set default starting values for nc 13 and nc14

for i=1:NumParticles
    if isempty(FitResults(i,1).Rate0)
        FitResults(i,1).Rate0=Rate014;
        FitResults(i,1).Elongation0=Elongation014;
        FitResults(i,1).TimeStart0=TimeStart014;
        FitResults(i,1).TimeEnd0=[];
        FitResults(i,1).RateOff0=[];
        FitResults(i,1).FrameFilter=[];
        FitResults(i,1).FitFrameRange=[];
        FitResults(i,1).Approved=0;
        FitResults(i,1).APID=0;
        FitResults(i,1).nc=[];
    elseif(isempty(FitResults(i,1).Elongation0))
        FitResults(i,1).Elongation0=Elongation014;
    end
end





         
%Go through each particles

FitFigure=figure;
% CurrentNC=12;
% i=min(find(sum(NParticlesAP)));
i=1;
cc=1;

while (cc~=13)  %13 is the ENTER button
    
    figure(FitFigure)
    clf
    
    if FitResults(i).Approved==-1
        set(gcf,'Color','r')
    elseif FitResults(i).Approved==1
        set(gcf,'Color','g')
    else
        set(gcf,'Color','default')
    end
    
    
    CurrentNC = CompiledParticles(i).nc;
    APbin = ceil(CompiledParticles(i).MeanAP/0.025);
    if APDivision(CurrentNC,APbin)
        if CurrentNC~=14
            FrameWindow=APDivision(CurrentNC,APbin):APDivision(CurrentNC+1,APbin);
        else
            FrameWindow=APDivision(CurrentNC,APbin):length(ElapsedTime);
        end
        NCFrameFilter = ismember(CompiledParticles(i).Frame,FrameWindow);

            %As an initial guess, use FrameFilter to determine the range of the
            %fit
            if isempty(FitResults(i).FitFrameRange)
                FitFrameRange = intersect(FrameWindow,CompiledParticles(i).Frame);
                if CurrentNC==14
                    FitFrameRange=FitFrameRange((ElapsedTime(FitFrameRange)-ElapsedTime(APDivision(CurrentNC,APbin)))<12);
                end
                FitResults(i).FitFrameRange=FitFrameRange;
            else
                FitFrameRange=FitResults(i).FitFrameRange;
            end
            
            FitFrameFilter=ismember(FitFrameRange,FrameWindow);
            %Extract the data for this range of frames
%             SlectedFrame = CompiledParticles(i).Frame(FitFrameFilter);
            FluoData=CompiledParticles(i).Fluo(NCFrameFilter);
            
            TimeData=ElapsedTime(CompiledParticles(i).Frame(NCFrameFilter));
            FluoDataForFit=FluoData(FitFrameFilter);
            FrameFilter = ismember(FrameWindow,CompiledParticles(i).Frame);

%             if CurrentNC<14
%                 %Do the fit
%                 x0=[FitResults(i).TimeStart0,...
%                     FitResults(i).TimeEnd0,...
%                     FitResults(i).Rate0,...
%                     FitResults(i).RateOff0,...
%                     FitResults(i).Elongation0];
%                 
%                 %set the boundary of parameter
%                 UpperBoundary = [30,50,1e4,1e4,3];
%                 LowerBoundary = [0,0,10,10,0.5];
% 
% 
%                 
%                 %Get rid of any NaN in the data
%                 NanFilter=~isnan(FluoData);
% 
%                 if ~isempty(TimeData(NanFilter))
% 
%                     [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
%                         lsqnonlin(@(x) lsqnonlinFitFluorescenceCurveV3(TimeData(NanFilter)-...
%                         ElapsedTime(FrameWindow(1)),...
%                         FluoData(NanFilter),GeneLength,...
%                         ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)),x),x0);
% 
%                     FitResults(i).TimeStart=xFit(1);
%                     FitResults(i).TimeEnd=xFit(2);
%                     FitResults(i).RateFit=xFit(3);
%                     FitResults(i).RateOffFit=xFit(4);
%                     FitResults(i).ElongationFit=xFit(5);
% 
%                     %Estimate an error bar out of the confidence intervals
%                     FitResults(i).CI=nlparci(xFit,residual,'jacobian',jacobian,'alpha',0.32); %modified from 0.68
% 
%                     FitResults(i).SDTimeStart=(FitResults(i).CI(1,2)-FitResults(i).CI(1,1))/2;
%                     FitResults(i).SDTimeEnd=(FitResults(i).CI(2,2)-FitResults(i).CI(2,1))/2;
%                     FitResults(i).SDRateFit=(FitResults(i).CI(3,2)-FitResults(i).CI(3,1))/2;
%                     FitResults(i).SDRateOffFit=(FitResults(i).CI(4,2)-FitResults(i).CI(4,1))/2;
%                     FitResults(i).SDElongationFit=(FitResults(i).CI(5,2)-FitResults(i).CI(5,1))/2;
% 
% 
%                     %Plot the results
%                     %Get the corresponding fitted curve
%                     [TimeFit,FluoFit]=FluorescenceCurveV3(ElapsedTime(FrameWindow(end))-...
%                         ElapsedTime(FrameWindow(1)),...
%                         xFit(1),xFit(2),xFit(3),xFit(4),xFit(5),GeneLength);
%                     
%                     %Plot the data could be used for the fit
%                     PlotHandle=plot(ElapsedTime(SlectedFrame)-ElapsedTime(FrameWindow(1)),FluoData(SlectedFrame),'.-k');
%                     hold on
%                     %Plot the data that was actually used for the fit
%                     PlotHandle(end+1)=plot(ElapsedTime(FitFrameRange)-ElapsedTime(FrameWindow(1)),...
%                         FluoData(ismember(SlectedFrame,FitFrameRange)),'or','MarkerFaceColor','r');
%                     
%                     %Plot the fit
%                     PlotHandle(end+1)=plot(TimeFit,FluoFit,'-r');
%                     hold off
%                     %StandardFigure(PlotHandle,gca)
%                     ylabel('Mean fluorescence nucleus')
%                     xlabel('Time into nc (min)')
%                     
%                     try
%                         ylim([0,max(MeanVectorAP(FrameWindow,i).*OnRatioAP(FrameWindow,i)/MaxOnRatio+...
%                             SDVectorAP(FrameWindow,i)./...
%                             sqrt(NParticlesAP(FrameWindow,i)).*OnRatioAP(FrameWindow,i)/MaxOnRatio)])
%                     catch
%                         display('Error in displaying the plot')
%                     end
% 
%                     legend(PlotHandle,['tON=',num2str(FitResults(i).TimeStart),' \pm ',num2str(FitResults(i).SDTimeStart)],...
%                         ['tOFF=',num2str(FitResults(i).TimeEnd),' \pm ',num2str(FitResults(i).SDTimeEnd)],...
%                         ['Rate=',num2str(FitResults(i).RateFit),' \pm ',num2str(FitResults(i).SDRateFit)],...
%                         ['RateOff=',num2str(FitResults(i).RateOffFit),' \pm ',num2str(FitResults(i).SDRateOffFit)],...
%                         ['Elongation=',num2str(FitResults(i).ElongationFit),' \pm ',num2str(FitResults(i).SDElongationFit)],...
%                         'Location','Best')
%                     
%   
%                 end
%             elseif CurrentNC==14
                
                %Do the fit
                x0=[FitResults(i).TimeStart0,FitResults(i).Rate0,FitResults(i).Elongation0];
                %set the boundary of parameter
                UpperBoundary14 = [30;1e4;3];
                LowerBoundary14 = [0;10;0.5];


                
                %Get rid of any NaN in the data
                NanFilter=~isnan(FluoDataForFit);

                if ~isempty(TimeData(NanFilter))

                    [xFit,resnorm,residual,exitflag,output,lambda,jacobian]=...
                        lsqnonlin(@(x) lsqnonlinFitSingleTrajectory(TimeData(NanFilter)-...
                        ElapsedTime(FrameWindow(1)),...
                        FluoDataForFit(NanFilter),GeneLength,...
                        ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)),x),x0);
                    
                    FitResults(i).APID = APbin;
                    FitResults(i).nc=CurrentNC;

                    FitResults(i).TimeStart=xFit(1);
                    FitResults(i).RateFit=xFit(2);
                    FitResults(i).ElongationFit=xFit(3);

                    %Estimate an error bar out of the confidence intervals
                    FitResults(i).CI=nlparci(xFit,residual,'jacobian',jacobian);

                    FitResults(i).SDTimeStart=(FitResults(i).CI(1,2)-FitResults(i).CI(1,1))/2;
                    FitResults(i).SDRateFit=(FitResults(i).CI(2,2)-FitResults(i).CI(2,1))/2;
                    FitResults(i).SDElongationFit=(FitResults(i).CI(3,2)-FitResults(i).CI(3,1))/2;


                    %Plot the results
                    %Get the corresponding fitted curve
                    [TimeFit,FluoFit]=FluorescenceCurveV0(ElapsedTime(FrameWindow(end))-...
                        ElapsedTime(FrameWindow(1)),...
                        xFit(1),1000,xFit(2),xFit(3),GeneLength);

                    %Plot the data that could be used for the fit
                    PlotHandle=plot(TimeData-ElapsedTime(FrameWindow(1)),...
                        FluoData,'o-k');
                    hold on
                    %Plot the data that was actually used for the fit
                    PlotHandle(end+1)=plot(ElapsedTime(FitFrameRange)-ElapsedTime(FrameWindow(1)),...
                        FluoData(FitFrameFilter),'or','MarkerFaceColor','r');
                    
                    %Plot the fit
                    PlotHandle(end+1)=plot(TimeFit,FluoFit,'-r');
                    hold off
                    ylabel('Mean fluorescence nucleus')
                    xlabel('Time into nc (min)')
                    
                    
                    try
                        ylim([0,max(MeanVectorAP(FrameWindow,i).*OnRatioAP(FrameWindow,i)/MaxOnRatio+...
                            SDVectorAP(FrameWindow,i)./...
                            sqrt(NParticlesAP(FrameWindow,i)).*OnRatioAP(FrameWindow,i)/MaxOnRatio)])
                    catch
                        display('Error in displaying the plot')
                    end

                    legend(PlotHandle,['tON=',num2str(FitResults(i).TimeStart),' \pm ',num2str(FitResults(i).SDTimeStart)],...
                        ['Rate=',num2str(FitResults(i).RateFit),' \pm ',num2str(FitResults(i).SDRateFit)],...
                        ['Elongation=',num2str(FitResults(i).ElongationFit),'\pm',num2str(FitResults(i).SDElongationFit)],...
                        'Location','Best')
                end

     end

    
    title([num2str(FitResults(i).APID),' AP, tStart0=',num2str(FitResults(i).TimeStart0),...
        ', tEnd0=',num2str(FitResults(i).TimeEnd0),', Rate=',num2str(FitResults(i).Rate0),...
        ', RateOff=',num2str(FitResults(i).RateOff0),', ElongationRate=',num2str(FitResults(i).Elongation0),...
        ', nc',num2str(CurrentNC)])
    
    
    %Set the limits on the x-axis
    if CurrentNC==14
%         xlim([0,ElapsedTime(end)])
        xlim([0,30])
%     else
%         xlim([0,ElapsedTime(eval(['nc',num2str(nc+1)]))-...
%             ElapsedTime(eval(['nc',num2str(nc)]))])
    end
    
    
    
    figure(FitFigure)
    ct=waitforbuttonpress;
    cc=get(FitFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    %Move between AP positions
    if (ct~=0)&(cc=='.')&(i<NumParticles)
        i=i+1;
    elseif (ct~=0)&(cc==',')&(i>1)
        i=i-1;
    
    %Approve, disapprove fit
    elseif (ct~=0)&(cc=='q')
        if FitResults(i).Approved==0
            FitResults(i).Approved=1;
        elseif FitResults(i).Approved==1
            FitResults(i).Approved=0;
        end

    
    %Disapprove, disapprove fit
    elseif (ct~=0)&(cc=='w')
        if FitResults(i).Approved==0
            FitResults(i).Approved=-1;
        elseif FitResults(i).Approved==-1
            FitResults(i).Approved=0;
        end
  
    
        
    %Move right range of fit
    elseif (ct~=0)&(cc=='k')&(length(FitResults(i).FitFrameRange)>2)
        FitResults(i).FitFrameRange=FitResults(i).FitFrameRange(1:end-1);
    elseif (ct~=0)&(cc=='l')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(i).FitFrameRange)))
            FilteredFramesTemp=FrameWindow(FrameFilter);
            FilterNotPrent = find(~ismember(FilteredFramesTemp,FitResults(i).FitFrameRange));
            SmallestLargerInx = min(find(FilteredFramesTemp(FilterNotPrent) > FitResults(i).FitFrameRange(end)));
            FitResults(i).FitFrameRange(end+1)= FilteredFramesTemp(FilterNotPrent(SmallestLargerInx));
%             FitResults(i).FitFrameRange(end+1)=...
%                 FilteredFramesTemp(min(find(~ismember(FilteredFramesTemp,FitResults(i).FitFrameRange))));
        end
    %Move left range of fit
    elseif (ct~=0)&(cc=='j')&(length(FitResults(i).FitFrameRange)>2)
        FitResults(i).FitFrameRange=FitResults(i).FitFrameRange(2:end);
    elseif (ct~=0)&(cc=='h')
        if ~isempty(find(~ismember(FrameWindow(FrameFilter),FitResults(i).FitFrameRange)))
            FilteredFramesTemp=FrameWindow(FrameFilter);
%             FilteredFramesTemp=FrameWindow(FrameFilter);
            FilterNotPrent = find(~ismember(FilteredFramesTemp,FitResults(i).FitFrameRange));
            LargestSmallerInx = max(find(FilteredFramesTemp(FilterNotPrent) < FitResults(i).FitFrameRange(1)));
%             FitResults(i).FitFrameRange(end+1)= FilteredFramesTemp(FilterNotPrent(LargestSmallerInx));
            FitResults(i).FitFrameRange=...
                [FilteredFramesTemp(FilterNotPrent(LargestSmallerInx)),FitResults(i).FitFrameRange];
%                 [FilteredFramesTemp(max(find(~ismember(FilteredFramesTemp,FitResults(i).FitFrameRange)))),...
%                 FitResults(i).FitFrameRange];
        end
    %Reset frame fit range
     elseif (ct~=0)&(cc=='r')   
        FitResults(i).FitFrameRange=FrameWindow(FrameFilter);

        
        
    %Change the initial parameters
    %TimeStart
    elseif (ct~=0)&(cc=='a')&((CurrentNC==14)|(CurrentNC~=14&FitResults(i).TimeStart0<FitResults(i).TimeEnd0))
        FitResults(i).TimeStart0=FitResults(i).TimeStart0+1;
    elseif (ct~=0)&(cc=='z')&(FitResults(i).TimeStart0>1)
        FitResults(i).TimeStart0=FitResults(i).TimeStart0-1;
    %TimeEnd
    elseif (ct~=0)&(cc=='s')&(FitResults(i).TimeEnd0<ElapsedTime(FrameWindow(end))-ElapsedTime(FrameWindow(1)))
        FitResults(i).TimeEnd0=FitResults(i).TimeEnd0+1;
    elseif (ct~=0)&(cc=='x')&(FitResults(i).TimeEnd0>FitResults(i).TimeStart0)
        FitResults(i).TimeEnd0=FitResults(i).TimeEnd0-1;
    %Rate, fine
    elseif (ct~=0)&(cc=='c')&(FitResults(i).Rate0>100)
        FitResults(i).Rate0=FitResults(i).Rate0-100;
    elseif (ct~=0)&(cc=='d')
        FitResults(i).Rate0=FitResults(i).Rate0+100;    
    %Rate, coarse
    elseif (ct~=0)&(cc=='C')&(FitResults(i).Rate0>100)
        FitResults(i).Rate0=FitResults(i).Rate0-500;
    elseif (ct~=0)&(cc=='D')
        FitResults(i).Rate0=FitResults(i).Rate0+500;   
    %Elongation Rate, coarse
    elseif (ct~=0)&(cc=='G')&(FitResults(i).Elongation0>0.5)
        FitResults(i).Elongation0=FitResults(i).Elongation0+0.25;
    elseif (ct~=0)&(cc=='B')
        FitResults(i).Elongation0=FitResults(i).Elongation0-0.25; 
    %Elongation rate, fine
    elseif (ct~=0)&(cc=='g')&(FitResults(i).Elongation0>0.5)
        FitResults(i).Elongation0=FitResults(i).Elongation0+0.05;
    elseif (ct~=0)&(cc=='b')
        FitResults(i).Elongation0=FitResults(i).Elongation0-0.05;    
    %RateOff, coarse
    elseif (ct~=0)&(cc=='V')&(FitResults(i).RateOff0<-100)
        FitResults(i).RateOff0=FitResults(i).RateOff0-500;
    elseif (ct~=0)&(cc=='F')
        FitResults(i).RateOff0=FitResults(i).RateOff0+500;
    %RateOff, fine
    elseif (ct~=0)&(cc=='v')&(FitResults(i).RateOff0<-100)
        FitResults(i).RateOff0=FitResults(i).RateOff0-100;
    elseif (ct~=0)&(cc=='f') 
        FitResults(i).RateOff0=FitResults(i).RateOff0+100;
        
        
    %Save
    elseif (ct~=0)&(cc=='e')
        save([DropboxFolder,filesep,Prefix,filesep,'SingleTrajectoryFits.mat'],...
        'FitResults')
    display('SingleTrajectoryFits.mat saved')
        
    %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard
    
    end
    
end

%Save the information
save([DropboxFolder,filesep,Prefix,filesep,'SingleTrajectoryFits.mat'],'FitResults')
display('FitResults.mat saved')        
        
close(FitFigure)
