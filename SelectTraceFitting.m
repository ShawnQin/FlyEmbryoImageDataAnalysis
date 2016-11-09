%This program select single particles traces
%Calculate the average fluorescence intensity that can be used to estimate loading
%rates and transcription start time

%revised on Sept 23,2016

%m,n    move to different nc
%,  .   move to different traces
%k, j   select the start and end frame
%s      save current data structure
%ENTER  exit


function SelectTraceFitting(varargin)
    
    %load the data
    Folder = '/Users/shan/Documents/MATLAB/LivemRNAFISH/Data/DynamicsResults';
    FolderTemp=uigetdir(Folder,'Choose folder with files to analyze');
    if(exist([FolderTemp,[filesep,'CompiledParticles.mat']])...
                && exist([FolderTemp,[filesep,'APDivision.mat']]))  
       load([FolderTemp,filesep,'CompiledParticles.mat']);
       load([FolderTemp,filesep,'APDivision.mat']);
    end
    
    %plot signle trace and selected certain region
    FrameNumThreshold = 15;  %only particle exits more than this frame is considered
    
    
    TraceFigure = figure;
    cc = 1;
    i0= 1;  %indice of particle trace
    CountTrace = 0;
    while(cc ~= 13)  %13 represent ENTER
        figure(TraceFigure)
        clf
%         while i0 <= length(CompiledParticles);
            ApprovedFrame = CompiledParticles(i0).Frame(CompiledParticles(i0).FrameApproved);
%             if(length(ApprovedFrame)>=FrameNumThreshold)
%                 plot(ApprovedFrame,CompiledParticles(i0).Fluo(CompiledParticles(i0).FrameApproved),'o-')
%             %select range
%             else
%                 i0 = i0+1;
%             end
%         end
        
        %plot and select
%         TimeWindow = ElapsedTime(ApprovedFrame);
        plot(ApprovedFrame,CompiledParticles(i0).Fluo(CompiledParticles(i0).FrameApproved),'o-')        
        ylabel('Fluorescence intensity')
        xlabel('Time into nc (min)')
        title([num2str(i0),'/',num2str(length(CompiledParticles))])
        
        
        
        figure(TraceFigure)
        ct=waitforbuttonpress;
        cc=get(TraceFigure,'currentcharacter');
        
        %move between different particle
        if (ct~=0)&&(cc=='.')&&(i0<length(CompiledParticles))
            i0 = i0 + 1;
        elseif(ct~=0)&&(cc==',')&&(i0>1)
            i0 = i0-1;
            
        %select certain frames
        elseif(ct ~=0)&&(cc=='j')
            [X1, Y1] = ginput(2);
            StartFrame = ceil(X1(1));
            EndFrame = floor(X1(2));
            [SelectedApprovedFrame,~, FrameInx]= intersect((StartFrame:EndFrame),ApprovedFrame);
%             FrameInx = ApprovedFrame(ApprovedFrame==SelectedApprovedFrame);
            CountTrace = CountTrace +1;
            SelectedTrace(CountTrace).Frame = SelectedApprovedFrame;
            SelectedTrace(CountTrace).Fluo = CompiledParticles(i0).Fluo(FrameInx);
            SelectedTrace(CountTrace).APbinID = ceil(CompiledParticles(i0).MeanAP/0.025);
            DivisionTimes = nnz(APDivision(:,SelectedTrace(CountTrace).APbinID));
            
            APIndex = SelectedTrace(CountTrace).APbinID;
            if (min(SelectedApprovedFrame)<APDivision(15-DivisionTimes,APIndex));
                    SelectedTrace(CountTrace).nc = 14-DivisionTimes;
            elseif(min(SelectedApprovedFrame)>=APDivision(14,APIndex));
                    SelectedTrace(CountTrace).nc = 14;
            else
                for j0 = 1:DivisionTimes-1;
                    if min(SelectedApprovedFrame) >= APDivision(14+j0-DivisionTimes,APIndex) && min(SelectedApprovedFrame)<APDivision(15+j0-DivisionTimes,APIndex)
                            SelectedTrace(CountTrace).nc = 14+j0-DivisionTimes;
                    end
                end
            end
            
            i0 = i0 + 1;
               
        %approve current trace without selection
        elseif(ct ~=0)&&(cc=='e')
            CountTrace = CountTrace+1;
            SelectedTrace(CountTrace).Frame = ApprovedFrame;
            SelectedTrace(CountTrace).Fluo = CompiledParticles(i0).Fluo;            
            SelectedTrace(CountTrace).APbinID = ceil(CompiledParticles(i0).MeanAP/0.025);
            SelectedTrace(CountTrace).nc = CompiledParticles(i0).nc;
            i0 = i0 + 1;
        %save
        elseif(ct ~=0)&&(cc=='s')
            save([FolderTemp,[filesep,'SelectedTrace.mat']],'SelectedTrace')
            display('SelectedTrace.mat saved')
        end
%         end
    end
end