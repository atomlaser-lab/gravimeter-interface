function varargout = testSequence(varargin)
    %% Initialize sequence - defaults should be handled here
    sq = initSequence;
    sq.find('3D Coils').set(0);
    sq.find('3D MOT Amp TTL').set(0);
    sq.find('25W Amp').set(0);
    sq.find('25W TTL').set(0);
    sq.find('25W Active').set(0);
    sq.find('50W TTL').set(0);
    sq.find('50W Pilot').set(1);
    sq.find('50W Amp').set(0);
    sq.delay(500e-3);
    sq.find('25W Active').set(1);
    sq.find('25W TTL').set(1);
    sq.find('50W TTL').set(1);
    T = 500e-3;
    t = 0:5e-3:T;
%     sq.find('25W Amp').after(t,sq.linramp(t,0,2));
    sq.find('50W Amp').after(t,sq.linramp(t,0,2));
    sq.delay(500e-3);
%     sq.find('25W Amp').set(sq.find('25W Amp').values(end) + 1);
    sq.find('50W Amp').set(sq.find('50W Amp').values(end) + 1);
%     sq.find('50W Amp').set(1);
    sq.delay(500e-3);
    sq.find('25W Amp').set(0);
    sq.find('50W Amp').set(0);
    sq.delay(100e-3);
    sq.find('25W TTL').set(0);
    sq.find('25W TTL').set(0);
    sq.find('25W Active').set(0);
    sq.find('50W TTL').set(0);
    sq.delay(100e-3);
    sq.find('25W TTL').set(0);
    sq.find('25W Active').set(0);
    sq.find('50W TTL').set(0);
    sq.find('50W Pilot').set(0);
    
    
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

