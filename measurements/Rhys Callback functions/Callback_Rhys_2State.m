function Callback_Rhys_2State(r)

% % Inputs
ClearImage = 0;
FigNum = 5;
TOF = 34e-3;

Title = '2 Pulse Raman: T = 1 ms, TOF = 16.5 ms, 1 AOM Power, 15 us Pulses, AOM detuning = 20.025 MHz';
Title = '1 Pulse Raman: TOF = 16.5 ms, 1 AOM Power, AOM detuning = 20.025 MHz';
Param = 0:2:100;
ParamName = 'Pulse Duration (us)';



if r.isInit()
    r.data.duration = Param;
    r.c.setup('var',r.data.duration);

elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.duration(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, T = %.0f us\n',r.c.now,r.c.total,...
        r.data.duration(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.5*rand);
    img = Abs_Analysis_DualState_RT('last');
    if ~img(1).raw.status.ok()
        
        Checks for an error in loading the files (caused by a missed
        image) and reruns the last sequence
        
        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end
    
    Store raw data
    
    r.data.files{i1,1} = img(1).raw.files;
    
    Get processed data
    
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');
    if numel(img) > 1
        r.data.N(i1,2) = r.data.N(i1,2) - r.data.N(i1,1);
        r.data.Nsum(i1,2) = r.data.Nsum(i1,2) - r.data.Nsum(i1,1);
    end
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));
    
    figure(FigNum);clf;
    subplot(1,2,1)
    scatter(r.data.duration(1:i1),r.data.N(1:i1,:),'filled'); %,'o'
    plot_format(ParamName,'Number','',12);

    title(' Raman frequency using fit over OD')
    grid on
    hold on;
    ylim([0,Inf]);
    
    subplot(1,2,2)
    scatter(r.data.duration(1:i1),r.data.Rsum(1:i1,:),'filled'); %,'o'
    hold off;
    plot_format(ParamName,'Population','',12);
    ylim([0,1])
    grid on;

    sgtitle(Title)
end