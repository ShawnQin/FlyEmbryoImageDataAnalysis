function [TimeRange,Fluorescence]=FluorescenceCurveV3(ncLength,TimeStart,TimeEnd,Rate,RateOff,Elongation,GeneLength)

%V3: Modified from V2, to have changing elongation rate
%V2: MOdified this to have an independent rate of turning off.

%This function outputs a standard fluorescent trace given a transcription
%rate, a start time, a delay time and an end time.




% TimeStart=5;
% TimeEnd=12;
% 
% Rate=4E3;     %Rate per minute
Delay = GeneLength/Elongation;
TimeRange=linspace(0,ncLength);

RateFrame=Rate*mean(diff(TimeRange));        %Rate per frame
RateOffFrame=RateOff*mean(diff(TimeRange));        %Rate per frame
Fluorescence=zeros(size(TimeRange));



for j=2:length(TimeRange)
    %Production part
    if (TimeRange(j)>TimeStart)&(TimeRange(j)<=(TimeStart+Delay))
        Fluorescence(j)=Fluorescence(j-1)+RateFrame;
    elseif (TimeRange(j)>(TimeStart+Delay))&(TimeRange(j)<=TimeEnd)
        Fluorescence(j)=Fluorescence(j-1);
  
    %Termination
    elseif (TimeRange(j)>TimeEnd)&(Fluorescence(j-1)>1E-10)%&(TimeRange(j)<=(-Rate*Delay/RateOff+TimeEnd))
        Fluorescence(j)=Fluorescence(j-1)+RateOffFrame;
    end
end
