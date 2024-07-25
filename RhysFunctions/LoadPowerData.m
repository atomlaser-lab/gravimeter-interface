
%% Inputs
Ch1FileLocation = "C:\Users\Apollo\Desktop\Temp\Ch2.csv";
Ch2FileLocation = "C:\Users\Apollo\Desktop\Temp\Ch1.csv";
FigNum = 20;
Title = '';
Title1 = 'Channel 1: Carrier';
Title2 = 'Channel 2: Sideband';

Row2Time = 30; %minutes
Row3Time = 60; %minutes

NormalisationPoint = 'mid';

MovingAverage = 20;

% YLIMITS = [0.9,1.3; 0.9, 1.15; 0.9, 1.15];
%% Load Ch1
% Load data
if exist('Ch1FileLocation','var') == 1
        warning('off','all')
    % % % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 4);

    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ["\t", " ", ",", ";", "_"];

    % Specify column names and types
    opts.VariableNames = ["VarName1", "VarName2", "PM", "VarName4"];
    opts.VariableTypes = ["datetime", "datetime", "categorical", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, "PM", "EmptyFieldRule", "auto");
    opts = setvaropts(opts, "VarName1", "InputFormat", "MM/dd/yyyy");
    opts = setvaropts(opts, "VarName2", "InputFormat", "HH:mm:ss");

    % Import the data
    tbl = readtable(Ch1FileLocation, opts);

    % % % Convert to output type
    VarName1 = tbl.VarName1;
    VarName2 = tbl.VarName2;
    VarName3 = tbl.PM;
    Ch1Amplitude = tbl.VarName4;

    % % % Clear temporary variables
    clear opts tbl
    warning('on','all')    
else
    Ch1Amplitude = [];
end

Ch1Amplitude = movmean(Ch1Amplitude,MovingAverage);

%% Load Ch2
if exist('Ch2FileLocation','var') == 1
    warning('off','all')    
    % Load data
    % % % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 4);

    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ["\t", " ", ",", ";", "_"];

    % Specify column names and types
    opts.VariableNames = ["VarName1", "VarName2", "PM", "VarName4"];
    opts.VariableTypes = ["datetime", "datetime", "categorical", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, "PM", "EmptyFieldRule", "auto");
    opts = setvaropts(opts, "VarName1", "InputFormat", "MM/dd/yyyy");
    opts = setvaropts(opts, "VarName2", "InputFormat", "HH:mm:ss");

    % Import the data
    tbl = readtable(Ch2FileLocation, opts);

    % % % Convert to output type
    VarName1 = tbl.VarName1;
    VarName2 = tbl.VarName2;
    VarName3 = tbl.PM;
    Ch2Amplitude = tbl.VarName4;

    % % % Clear temporary variables
    clear opts tbl
    warning('on','all')
else
    Ch2Amplitude = [];
end

Ch2Amplitude = movmean(Ch2Amplitude,MovingAverage);


%% Make nan values
if exist('Ch2FileLocation','var') == 0
    Ch2Amplitude =  nan(size(Ch1Amplitude));
elseif exist('Ch1FileLocation','var') == 0
    Ch1Amplitude = nan(size(Ch2Amplitude));
end

clear Ch1FileLocation Ch2FileLocation
%% Convert time data to seconds
Var_24Hour = VarName2 + hours(VarName3 == "PM" & VarName2.Hour~=12)*12;
HourChangeIndex = find(abs(diff(Var_24Hour.Hour)) >= 1 == 1);
NewHours = zeros(size(Var_24Hour.Hour));

if numel(HourChangeIndex) ~= 1
    for ii = 1:numel(HourChangeIndex)
        if ii == 1
            NewHours(1:HourChangeIndex(ii)) = Var_24Hour.Hour(1);
        elseif ii == numel(HourChangeIndex)
            NewHours(HourChangeIndex(end-1)+1:HourChangeIndex(end)) = NewHours(HourChangeIndex(end-1)) + 1;
            NewHours(HourChangeIndex(end)+1 : end) = NewHours(HourChangeIndex(end)) + 1;
        else
            NewHours(HourChangeIndex(ii-1)+1:HourChangeIndex(ii)) = NewHours(HourChangeIndex(ii-1)) + 1;
        end
    end
else
    NewHours(1:HourChangeIndex) = Var_24Hour.Hour(1);
    NewHours(HourChangeIndex+1:end) = Var_24Hour.Hour(1) + 1;
end
NewHours = NewHours-NewHours(1);

TotalTime_seconds = NewHours*60*60 + Var_24Hour.Minute*60 + Var_24Hour.Second;
TotalTime_seconds = TotalTime_seconds - TotalTime_seconds(1);

%% Plot


Row2Index = find(TotalTime_seconds - Row2Time*60 > 0,1);
Row3Index = find(TotalTime_seconds - Row3Time*60 > 0,1);
if isstr(NormalisationPoint) == 1
    if strcmpi(NormalisationPoint,'end') == 1
        NormalisationPoint = length(Ch1Amplitude);
    elseif strcmpi(NormalisationPoint,'mid') == 1
        NormalisationPoint = round(length(Ch1Amplitude)/2,0);        
    end
end


figure(FigNum);clf
sgtitle(Title)
subplot(3,2,1)
scatter(TotalTime_seconds/60/60,Ch1Amplitude/(Ch1Amplitude(NormalisationPoint)),5,'k')
hold on
plot(TotalTime_seconds/60/60,Ch1Amplitude/(Ch1Amplitude(NormalisationPoint)),'b','LineWidth',1)
xlim([0,TotalTime_seconds(end)/60/60])
% plot(TotalTime_seconds/60/60,Ch1Amplitude,'b')
title(Title1)
xlabel('time (hours)')
ylabel('Amp/Amp_{0}')
if exist('YLIMITS','var') == 1
    ylim(YLIMITS(1,:))
end

subplot(3,2,2)
scatter(TotalTime_seconds/60/60,Ch2Amplitude/(Ch2Amplitude(NormalisationPoint)),5,'k')
hold on
plot(TotalTime_seconds/60/60,Ch2Amplitude/(Ch2Amplitude(NormalisationPoint)),'r')
xlim([0,TotalTime_seconds(end)/60/60])

title(Title2)
xlabel('time (hours)')
ylabel('Amp/Amp_{0}')
if exist('YLIMITS','var') == 1
    ylim(YLIMITS(1,:))
end

subplot(3,2,3)
scatter(TotalTime_seconds(1:Row2Index)/60,Ch1Amplitude(1:Row2Index)/(Ch1Amplitude(NormalisationPoint)),5,'k')
hold on
plot(TotalTime_seconds(1:Row2Index)/60,Ch1Amplitude(1:Row2Index)/(Ch1Amplitude(NormalisationPoint)),'b')
xlabel('time (min)')
ylabel('Amp/Amp_{0}')
xlim([0,Row2Time])
if exist('YLIMITS','var') == 1
    ylim(YLIMITS(2,:))
end

subplot(3,2,4)
scatter(TotalTime_seconds(1:Row2Index)/60,Ch2Amplitude(1:Row2Index)/(Ch2Amplitude(NormalisationPoint)),5,'k')
hold on
plot(TotalTime_seconds(1:Row2Index)/60,Ch2Amplitude(1:Row2Index)/(Ch2Amplitude(NormalisationPoint)),'r')
xlabel('time (min)')
ylabel('Amp/Amp_{0}')
xlim([0,Row2Time])
if exist('YLIMITS','var') == 1
    ylim(YLIMITS(2,:))
end


subplot(3,2,5)
scatter(TotalTime_seconds(1:Row3Index)/60,Ch1Amplitude(1:Row3Index)/(Ch1Amplitude(NormalisationPoint)),5,'k')
hold on
plot(TotalTime_seconds(1:Row3Index)/60,Ch1Amplitude(1:Row3Index)/(Ch1Amplitude(NormalisationPoint)),'b')
xlabel('time (min)')
ylabel('Amp/Amp_{0}')
xlim([0,Row3Time])
if exist('YLIMITS','var') == 1
    ylim(YLIMITS(3,:))
end

subplot(3,2,6)
scatter(TotalTime_seconds(1:Row3Index)/60,Ch2Amplitude(1:Row3Index)/(Ch2Amplitude(NormalisationPoint)),5,'k')
hold on
plot(TotalTime_seconds(1:Row3Index)/60,Ch2Amplitude(1:Row3Index)/(Ch2Amplitude(NormalisationPoint)),'r')
xlabel('time (min)')
ylabel('Amp/Amp_{0}')
xlim([0,Row3Time])
if exist('YLIMITS','var') == 1
    ylim(YLIMITS(3,:))
end

PowerRatio = Ch1Amplitude./(Ch2Amplitude);
figure(FigNum+1);clf
subplot(3,1,1)
scatter(TotalTime_seconds/60/60,PowerRatio/PowerRatio(NormalisationPoint),5,'k')
hold on
plot(TotalTime_seconds/60/60,PowerRatio/PowerRatio(NormalisationPoint),'r')
xlabel('time (hours)')
ylabel('I2/I1')
xlim([0 TotalTime_seconds(end)/60/60])
subplot(3,1,2)
scatter(TotalTime_seconds(1:Row2Index)/60,PowerRatio(1:Row2Index)/PowerRatio(NormalisationPoint),5,'k')
hold on
plot(TotalTime_seconds(1:Row2Index)/60,PowerRatio(1:Row2Index)/PowerRatio(NormalisationPoint),'r')
xlabel('time (min)')
ylabel('I2/I1')
xlim([0 TotalTime_seconds(Row2Index)/60])
subplot(3,1,3)
scatter(TotalTime_seconds(1:Row3Index)/60,PowerRatio(1:Row3Index)/PowerRatio(NormalisationPoint),5,'k')
hold on
plot(TotalTime_seconds(1:Row3Index)/60,PowerRatio(1:Row3Index)/PowerRatio(NormalisationPoint),'r')
xlim([0 TotalTime_seconds(Row3Index)/60])
xlabel('time (min)')
ylabel('I2/I1')

clear PowerRatio VarName1 VarName2 VarName3 Var_24Hour YLIMITS Title2 Title1 Title Row3Time Row3Index Row2Time Row2Index NewHours ii HourChangeIndex FigNum