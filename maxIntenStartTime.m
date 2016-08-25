function [MeanMaxIntensityAP,MeanStartTimeAP,SelectedTrace,MeanmRNAAP] ...
    = maxIntenStartTime(varargin)
if length(varargin)==1
    FolderTemp=varargin{1};
    PlotFlag = 0;
elseif(length(varargin>1))
    PlotFlag = varargin{2}; %determine if plot the figure
else
    PlotFlag = 0;
end
    
    
%load the data
load([FolderTemp,filesep,'CompiledParticles.mat']);
load([FolderTemp,filesep,'APDivision.mat']);

%maximum intensity of each particle
% MaxIntensityParticle = zeros(length(CompiledParticles),3); %with the second column as the nc index and third column as AP index
MaxIntensityParticle = [];
mRNAParticle = [];
ParticleStartTimeDuration = [];
ParticleSelected = [];
SelectedTrace = [];
MaxFrame = size(AllTracesVector,1);
FrameNumThreshold = 15;  %only particle exits more than this frame is considered
for i0 = 1:length(CompiledParticles);
%     MaxIntensityParticle = nanmax(AllTracesVector);
    ApprovedFrame = CompiledParticles(i0).Frame(CompiledParticles(i0).FrameApproved);
    if(length(ApprovedFrame)>=FrameNumThreshold)
        MaxIntensityParticle = [MaxIntensityParticle;max(CompiledParticles(i0).Fluo),...
            CompiledParticles(i0).nc,find(APFilter(i0,:)==1)];
        ParticleSelected = [ParticleSelected,i0];
        SelectedTrace = [SelectedTrace,nan(MaxFrame,1)];
        
        %some particles across more than one nc, should consider this situation
        %somehow the direct read of nc is not acurate, here i would like to
        %recalculate it base on the APDivision and ApprovedFrame
        NCnumber  = sum(APDivision(:,find(APFilter(i0,:)==1))~=0); %measured nc number
        NCFrame = APDivision(APDivision(:,find(APFilter(i0,:)==1))~=0,find(APFilter(i0,:)==1));
        AllNCInterval = cell(NCnumber+1,1);
        AllNCInterval{1} =  1:NCFrame(1);
        AllNCInterval{NCnumber + 1} = NCFrame(end):size(AllTracesVector,1);
        for k0 = 1:NCnumber-1;
            AllNCInterval{k0+1} = NCFrame(k0):NCFrame(k0+1);
        end
        
        AllCoverRatio = zeros(NCnumber,1);
       
        for k0 = 1:NCnumber+1;
            AllCoverRatio(k0) = sum(ismember(ApprovedFrame,AllNCInterval{k0}))/length(AllNCInterval{k0});
        end
        [~,WhichNC] = nanmax(AllCoverRatio);
        ncFrame = ApprovedFrame(ismember(ApprovedFrame,AllNCInterval{WhichNC}));
        MaxIntensityParticle(end,2) = 13-NCnumber + WhichNC;  %modify the nc index of this particle
%         if(ApprovedFrame(end)<=NCFrame(1))
%             ncDuration = 1:NCFrame(1);
%             MaxIntensityParticle(end,2) = 14-NCnumber; %modify the nc
%         elseif(ApprovedFrame(1)>=NCFrame(end))
%             ncDuration = NCFrame(end):max(ApprovedFrame);
%             MaxIntensityParticle(end,2) = 14;
%         else
%             for j0 = 1:NCnumber-1;
%                 if(ApprovedFrame(1)>=NCFrame(j0) && ApprovedFrame(1)<NCFrame(j0+1))
%                    ncDuration  = NCFrame(j0):NCFrame(j0+1);
%                    MaxIntensityParticle(end,2) = 14- NCnumber + i0;
%                    continue
%                 end
%             end
%         end
            
            
%         if CompiledParticles(i0).nc ~=14
%             if(APDivision(CompiledParticles(i0).nc,find(APFilter(i0,:)==1))==0)
%                 if(min(ApprovedFrame)>=APDivision(CompiledParticles(i0).nc+1,find(APFilter(i0,:)==1)))
%                     ncDuration = APDivision(CompiledParticles(i0).nc+1,find(APFilter(i0,:)==1)):APDivision(CompiledParticles(i0).nc+2,find(APFilter(i0,:)==1));
%                 else
%                     ncDuration = min(ApprovedFrame):APDivision(CompiledParticles(i0).nc+1,find(APFilter(i0,:)==1));
%                 end
%             else               
%                 ncDuration = APDivision(CompiledParticles(i0).nc,find(APFilter(i0,:)==1)):APDivision(CompiledParticles(i0).nc+1,find(APFilter(i0,:)==1));
%             end
%         else
%             ncDuration = APDivision(CompiledParticles(i0).nc,find(APFilter(i0,:)==1)):max(ApprovedFrame);
%         end
%         ncFrame = ApprovedFrame(ismember(ApprovedFrame,ncDuration));
%         startInx = min(ApprovedFrame(ismember(ApprovedFrame,ncDuration)));
%         terminationInx = max(ApprovedFrame(ismember(ApprovedFrame,ncDuration)));
        if(ncFrame(1)~=AllNCInterval{WhichNC}(1))
            ParticleStartTimeDuration = [ParticleStartTimeDuration;ElapsedTime(ncFrame(1))- ElapsedTime(AllNCInterval{WhichNC}(1)),...
            MaxIntensityParticle(end,2)];
        else
            %particles that across different nc
            ParticleStartTimeDuration = [ParticleStartTimeDuration;nan,MaxIntensityParticle(end,2)];
        end
        if(length(ncFrame)>1)
            mRNAParticle = [mRNAParticle;trapz(ElapsedTime(ncFrame),AllTracesVector(ncFrame,i0)),MaxIntensityParticle(end,2)];
        else
            mRNAParticle = [mRNAParticle;nan,MaxIntensityParticle(end,2)];
        end
        SelectedTrace(ncFrame,end) = AllTracesVector(ncFrame,i0);
%         MaxIntensityParticle(i0,2) = CompiledParticles(i0).nc;
%         MaxIntensityParticle(i0,3) = find(APFilter(i0,:)==1);
    end
end

%maximum intensity along AP
MeanMaxIntensityAP = nan(4,41);
StdMaxIntensityAP = nan(4,41);
MeanStartTimeAP = nan(4,41);
StdStartTimeAP = nan(4,41);
MeanmRNAAP = nan(4,41);
StdmRNAAP = nan(4,41);

APInx = unique(MaxIntensityParticle(:,3));
for i0 = 1:4;
    for j0 = 1:length(unique(MaxIntensityParticle(:,3)))
        SelectInx = MaxIntensityParticle(:,2)==(i0+10) & MaxIntensityParticle(:,3)==APInx(j0);
        MeanMaxIntensityAP(i0,APInx(j0)) = mean(MaxIntensityParticle(SelectInx,1));
        StdMaxIntensityAP(i0,APInx(j0)) = std(MaxIntensityParticle(SelectInx,1));
        
        MeanStartTimeAP(i0,APInx(j0)) = nanmean(ParticleStartTimeDuration(SelectInx,1));
        StdStartTimeAP(i0,APInx(j0)) = nanstd(ParticleStartTimeDuration(SelectInx,1));
        
        MeanmRNAAP(i0,APInx(j0)) = nanmean(mRNAParticle(SelectInx,1));
        StdmRNAAP(i0,APInx(j0)) = nanstd(mRNAParticle(SelectInx,1));
    end
end


if PlotFlag
%plot all selected single particles
figure(1)
plot(ElapsedTime,SelectedTrace)
xlabel('Time(min)','FontSize',24,'FontWeight','Bold')
ylabel('Fluorescence intensity(a.u)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')

%plot the mean profile compared without selection
figure(2)
hold on
plot(ElapsedTime,nanmean(SelectedTrace,2))
plot(ElapsedTime,nanmean(AllTracesVector,2))
hold off
xlabel('Time(min)','FontSize',24,'FontWeight','Bold')
ylabel('Fluorescence intensity(a.u)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')


%plot the average transcription start time along AP
figure(3)
hold on
for i0 = 1:4
    errorbar((0:0.025:1),MeanStartTimeAP(i0,:),StdStartTimeAP(i0,:))
end
hold off
xlabel('AP position','FontSize',24,'FontWeight','Bold')
ylabel('Transcription start time(min)','FontSize',24,'FontWeight','Bold')
set(gca,'FontSize',24,'FontWeight','Bold')
end

end