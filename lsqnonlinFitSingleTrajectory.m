function chi2=lsqnonlinFitSingleTrajectory(TimeData,FluoData,GeneLength,ncLength,x0)

%V2 this version of also fit the elongation rate
%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin. This is for cycle 14, where we don't fit for the signal going
%down.

%Starting conditions
TimeStart0=x0(1);
Rate0=x0(2);
Elongation0 = x0(3);

Delay = GeneLength/Elongation0;

for i=1:length(TimeData)
    if TimeData(i)<=TimeStart0
        FluoPrediction(i)=0;
    elseif (TimeData(i)>TimeStart0)&(TimeData(i)<=TimeStart0+Delay)
        FluoPrediction(i)=Rate0*(TimeData(i)-TimeStart0);
    elseif (TimeData(i)>TimeStart0+Delay)
        FluoPrediction(i)=Rate0*Delay;
    end
end



%[TimeRange,Fluorescence]=FluorescenceCurve(ncLength,TimeStart0,TimeEnd0,Rate0,Delay);


%Interpolate the values of the generated curve for the actual time points
%we have
%InterpData=pchip(TimeRange,Fluorescence,TimeData);

%NanFilter=~isnan(FluoData);
%chi2=(FluoData(NanFilter)-InterpData(NanFilter)').^2;
try
    chi2=(FluoData-FluoPrediction).^2;
catch
    1+1;
end

