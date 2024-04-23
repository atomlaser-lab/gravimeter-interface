function varargout = testDDSTables(varargin)
    %% Initialize sequence - defaults should be handled here
    sq = initSequence;

%     sq.dds(1).power_conversion_method = DDSChannel.POWER_CONVERSION_DBM_INTERP;
%     sq.dds(2).power_conversion_method = DDSChannel.POWER_CONVERSION_DBM_INTERP;
%     calibData = load('RamanAOMData_11042024');
%     sq.dds(1).calibrationData = calibData.data_ch1;
%     sq.dds(2).calibrationData = calibData.data_ch2;

    sq.dds(1).set(20,0,0);
    sq.dds(2).set(20,0,0);
    sq.delay(500e-3);
    timeAtDrop = sq.time;
    

    %% Interferometry
    % Issue falling-edge trigger for MOGLabs DDS box
    sq.find('raman dds trig').before(10e-3,1);
    sq.find('raman dds trig').after(10e-3,0); %MOGLabs DDS triggers on falling edge
    sq.find('raman dds trig').after(50e-3,1);
    
    % Create a sequence of Bragg pulses. The property ddsTrigDelay is used
    % in compiling the DDS instructions and making sure that they start at
    % the correct time.
    sq.ddsTrigDelay = timeAtDrop;   
    t = 0:1:20;
    P = 0.1*ones(size(t));
    freq = 110*ones(size(t));
    ph = 0*ones(size(t));
    
    sq.dds(1).after(t*1e-6 + 5e-6,freq,P,ph);
    sq.dds(2).after(t*1e-6 + 5e-6,freq,P,ph);

    sq.anchor(sq.latest);
    sq.delay(1e-3);
    sq.dds(1).set(20,0,0);
    sq.dds(2).set(20,0,0);
    sq.delay(1e-3);
    
    
    %% Automatic start
    %If no output argument is requested, then compile and run the above
    %sequence
    if nargout == 0
        r = RemoteControl;
        r.upload(sq.compile);
        r.run;
    else
        varargout{1} = sq;
    end

end

