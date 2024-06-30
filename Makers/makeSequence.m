function varargout = makeSequence(varargin)

%% New parse inputs
% %% comment out below parse inputs for this old script to work
allVarNames = evalin('base', 'who');
seqOptVarIndices = cellfun(@(varName) isa(evalin('base', varName), 'SequenceOptions'), allVarNames);
seqOptVarNames = allVarNames(seqOptVarIndices);

if isempty(seqOptVarNames)
    if  nargin ==1 && ~exist(varargin{1},'var')
        error('define the sequence option variable name first');
    end
    % Create a new instance of SequenceOptions if it does not exist
    opt = SequenceOptions;
    cprintf('Keywords','Initialising Sequence Options\n')
    assignin('base', 'opt', opt);
else
    if numel(seqOptVarNames) > 1
        error('Multiple SequenceOptions found in the workspace');
    end
    % Use the existing instance of SequenceOptions
    opt = evalin('base', seqOptVarNames{1});
end

%% Parse input arguments
% opt = SequenceOptions;
% 
% if nargin == 1
%     %
%     % If first argument is of type GravimeterOptions, use that
%     %
%     if ~isa(varargin{1},'SequenceOptions')
%         error('If using one argument it must be of type SequenceOptions');
%     end
%     opt = opt.replace(varargin{1}); 
% elseif mod(nargin,2) == 0
%     opt.set(varargin{:});
% elseif mod(nargin - 1,2) == 0 && isa(varargin{1},'SequenceOptions')
%     opt = opt.replace(varargin{1});
%     opt.set(varargin{2:end});    
% else 
%     error('You must supply either a SequenceOptions object, a set of name/value pairs for options, or a SequenceOptions object followed by name/value pairs');
% end

%% Create conversion object
convert = RunConversions;
imageVoltage = convert.imaging(opt.detuning - 1.8);
%% Create a BEC
sq = makeBEC(opt);

%% Measure trap freq.
% sq.find('50W Amp').set(0).after(0.25e-3,sq.find('50W Amp').values(end-1));
% % sq.find('25W Amp').set(convert.dipole25(opt.dipoles + 0.0));
% sq.delay(0.5e-3 + opt.params(1));

%% Drive width oscillations
% f = 71;
% Nosc = 5;
% T = Nosc/f;
% t = linspace(0,T,20*Nosc);
% sq.find('50W amp').after(t,convert.dipole50(opt.dipoles + 0.15*sin(2*pi*f*t)));
% sq.delay(max(t) + opt.params(1));

%% Drop atoms
% sq.find('mot coil ttl').set(0);
% sq.find('3D Coils').set(convert.mot_coil(0));
% sq.delay(35e-3 - opt.tof);

timeAtDrop = sq.time; %Store the time when the atoms are dropped for later
sq.anchor(timeAtDrop);
sq.find('3D mot amp ttl').set(0);
% sq.find('bias e/w').before(200e-3,0);
% sq.find('bias n/s').before(200e-3,0);
% sq.find('bias u/d').before(200e-3,0);
sq.find('mw amp ttl').set(0);
sq.find('mot coil ttl').set(0);
sq.find('3D Coils').set(convert.mot_coil(0));
sq.find('25w ttl').set(0);
sq.find('50w ttl').set(0);
sq.find('50w amp').set(convert.dipole50(0));
sq.find('25w amp').set(convert.dipole25(0)); 

%% Interferometry
enableDDS = 1;      %Enable DDS and DDS trigger
enableBragg = 1;    %Enable Bragg diffraction
enableRaman = 0;    %Enable Raman transition
enableGrad = 0;     %Enable gradiometry
if enableDDS
    % 
    % Issue falling-edge trigger for MOGLabs DDS box when DDS is
    % enabled
    %
    sq.find('dds trig').before(10e-3,1);
    sq.find('dds trig').after(10e-3,0); %MOGLabs DDS triggers on falling edge
    sq.find('dds trig').after(10e-3,1);
    sq.ddsTrigDelay = timeAtDrop; 
end

if enableDDS && enableBragg
    %
    % Create a sequence of Bragg pulses. The property ddsTrigDelay is used
    % in compiling the DDS instructions and making sure that they start at
    % the correct time.
    %
%         sq.find('liquid crystal bragg').at(sq.find('Liquid crystal Bragg').last,-2); %Use this to get a little more light on the phase lock PD
    braggOrder = 1;
    k = 2*pi*384229441689483/const.c;  %Frequency of Rb-85 F=3 -> F'=4 transition
    vrel = abs(2*const.hbar*k/const.mRb);
    dv = 700e-6/216.5e-3;
    chirp = opt.bragg.chirp;
    T = opt.bragg.T;
    Tasym = opt.bragg.Tasym;
    %
    % Calculate seperation time - this depends only on the asymmetry,
    % bragg order, and time of flight.
    %
    if Tasym == 0
        dsep = 2*dv*opt.tof;
        Tsep = dsep/vrel;
    else  
        Tsep = 1*abs(const.mRb*pi*opt.tof/(4*braggOrder*k^2*const.hbar*Tasym));
    end
    % Override the calculated separation time if specified by user
    if ~isempty(opt.bragg.Tsep)
        Tsep = opt.bragg.Tsep;
    end
    %
    % Calculate when the initial pulse should arrive
    %
    if enableGrad
        if isempty(opt.bragg.t0)
            t0 = 30e-3;
        else
            t0 = opt.bragg.t0;
        end

        if isempty(opt.bragg.ti)
            ti = opt.tof - t0 - Tasym - Tsep - 2*T;
        else
            ti = opt.bragg.ti;
        end
    else
        if isempty(opt.bragg.t0)
            t0 = opt.tof - 2*T - Tsep - Tasym;
        else
            t0 = opt.bragg.t0;
            Tsep = opt.tof - 2*T - t0;
        end
    end

%         if numel(t0) > 1
%             error('Unable to determine t0!');
%         elseif t0 < 30e-3
%             warning('Initial Bragg pulse occurs at %.1f ms and will be clamped to 30 ms!',t0*1e3);
%         end
%         t0 = max(t0,30e-3);
    if enableGrad
        fprintf(1,'t0 = %0.3f ms, ti = %0.3f ms, Tsep = %0.3f ms\n',t0*1e3,ti*1e3,Tsep*1e3);
    else
        fprintf(1,'t0 = %0.3f ms, Tsep = %0.3f ms\n',t0*1e3,Tsep*1e3);
    end

    if enableGrad
        %
        % Initial velocity selection
        %
        makeBraggSequence(sq.dds,'k',k,'dt',1e-6,'t0',t0,'T',T,...
            'width',60e-6,'Tasym',0,'phase',[0,0,0],'chirp',chirp,...
            'power',[0.6,0,0],'order',-4);
        %
        % Splitting of the cloud
        %
        sq.dds.anchor(timeAtDrop);
        makeBraggSequence(sq.dds,'k',k,'dt',1e-6,'t0',20e-3 + t0,'T',T,...
            'width',opt.bragg.width,'Tasym',Tasym,'phase',0,'chirp',chirp,...
            'power',0.35.*0,'order',braggOrder,'start_order',-4);
        %
        % Interferometry
        %
        sq.dds.anchor(timeAtDrop);
        makeBraggSequence(sq.dds,'k',k,'dt',1e-6,'t0',20e-3 + t0 + ti,'T',T,...
            'width',opt.bragg.width,'Tasym',Tasym,'phase',[0,0,opt.bragg.phase],'chirp',chirp,...
            'power',opt.bragg.power.*[1,1,1]*0,'order',braggOrder,'start_order',-4);
    else
        makeBraggSequence(sq.dds,'k',k,'dt',1e-6,'t0',t0,'T',T,...
            'width',opt.bragg.width,'Tasym',Tasym,'phase',[0,0,opt.bragg.phase],'chirp',chirp,...
            'power',opt.bragg.power*[1,2,1],'order',braggOrder);

%             makeCompositePulse(sq.dds,'k',k,'dt',1e-6,'t0',t0,'pulse separation',100e-6,...
%                 'width',30e-6,'phase',[0,180,0],'chirp',chirp,...
%                 'power',opt.bragg_power*[90/180,0/180,0/180],'order',braggOrder);

%             [amp,ph,freq,flags,t] = makeBraggSequence_pl('t0',t0,'T',T,'fwhm',40e-6,...
%                 'power',opt.bragg_power*[1,0,0],'phase',2*[0,0,opt.final_phase]*pi/180,...
%                 'usehold',0,'holdfreq',2,'holdamp',0.04);
%             pl = PhaseLock('192.168.1.38');
%             pl.setDefaults;
%             pl.amp(1).set(0);
%             pl.amp(2).set(0);
%             pl.useManual.set(0);
%             pl.disableExtTrig.set(0);
%             pl.shift.set(3);
%             pl.upload;
%             pl.uploadTiming(t,ph,amp,freq,flags);

    end


end

if enableDDS && enableRaman
    %
    % Start by re-anchoring the internal pointer for the DDS channels at the
    % drop time.  I do this to make referencing the time at which the Raman
    % pulse occurs easier to calculate
    %
    sq.dds.anchor(timeAtDrop);
    sq.dds(1).after(1e-4,DDSChannel.DEFAULT_FREQ - opt.raman.df/4,0,0);
    sq.dds(2).after(1e-4,DDSChannel.DEFAULT_FREQ + opt.raman.df/4,0,0);
    sq.delay(50e-3);
    T = opt.raman.width;
    dt = 1e-6;
    t = 0:dt:T;
    sq.dds(1).after(t,DDSChannel.DEFAULT_FREQ - opt.raman.df/4,opt.raman.power,0);
    sq.dds(2).after(t,DDSChannel.DEFAULT_FREQ + opt.raman.df/4,opt.raman.power,0);
    sq.delay(T);
    sq.dds(1).after(t,DDSChannel.DEFAULT_FREQ,0,0);
    sq.dds(2).after(t,DDSChannel.DEFAULT_FREQ,0,0);
end

if any(opt.mw.enable)
    %
    % Apply a pair of microwave pulses to effect the transfers
    % |F=1,m=-1> -> |F=2,m=0> -> |F=1,m=0>.  The first pulse is applied
    % 10 ms after the atoms are dropped to minimize any possible
    % state-changing collisions.  The "R&S list step trig" skips to the
    % next frequency on the rising edge and resets the list on the
    % falling edge
    %
    if opt.mw.enable(1)
        sq.anchor(timeAtDrop);
        sq.find('bias e/w').set(10);
        sq.find('R&S list step trig').set(1);
        sq.delay(20e-3);
        sq.find('state prep ttl').set(1);
        sq.delay(372e-6);
        sq.find('state prep ttl').set(0);
    end

    if opt.mw.enable(2) && ~opt.mw.analyze(1)
        sq.find('Repump Amp TTL').set(1).after(1e-3,0);
        sq.find('Liquid Crystal Repump').set(-2.22).after(1e-3,7);
        sq.find('repump freq').set(4.3); 

        sq.find('R&S list step trig').set(0);
        sq.delay(20e-3);
        sq.find('state prep ttl').set(1);
        sq.delay(250e-6);
        sq.find('state prep ttl').set(0);

        sq.find('R&S list step trig').set(1);
        sq.find('bias e/w').set(0);
        sq.find('3D MOT Amp TTL').set(1).after(100e-6,0);
    end
else
    sq.find('bias e/w').at(timeAtDrop,0);
end

if opt.mw.enable_sg || opt.mw.analyze(1)
    %
    % Apply a Stern-Gerlach pulse to separate states based on magnetic
    % moment.  A ramp is used to ensure that the magnetic states
    % adiabatically follow the magnetic field
    %     
    sq.anchor(timeAtDrop + 50e-3);
%     sq.delay(30e-3);
%     sq.waitFromLatest(5e-3);
    Tsg = 5e-3;
    sq.find('mot coil ttl').set(1);
    t = linspace(0,Tsg,20);
    sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,0,1.5)));
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
    sq.delay(2*Tsg);
    sq.find('mot coil ttl').set(0);
    sq.find('3d coils').set(convert.mot_coil(0));
end

%% Imaging stage
%
% Image the atoms.  Reset the pointer for the whole sequence to when
% the atoms are dropped from the trap.  This means that the
% time-of-flight (tof) used in makeImagingSequence is now the delay
% from the time at which the atoms are dropped to when the first
% imaging pulse occurs
%
% sq.anchor(timeAtDrop);
% sq.camDelay = timeAtDrop - 2;   %Set camera acquisition delay to be 2 s less than when image is taken
% if strcmpi(opt.imaging_type,'drop 1') || strcmpi(opt.imaging_type,'drop 2')
%     makeImagingSequence(sq,'type',opt.imaging_type,'tof',opt.tof,...
%         'repump Time',100e-6,'pulse Delay',10e-6,'pulse time',[],...
%         'imaging freq',imageVoltage,'repump delay',10e-6,'repump freq',4.3,...
%         'manifold',1,'includeDarkImage',true,'cycle time',150e-3);
% elseif strcmpi(opt.imaging_type,'drop 3') || strcmpi(opt.imaging_type,'drop 4')
%     makeFMISequence(sq,'tof',opt.tof,'offset',30e-3,'duration',100e-3,...
%         'imaging freq',imageVoltage,'manifold',1);
% end

%% Automatic save of run
%
% This automatically creates a record of this file when the sequence is
% created
%
fpathfull = [mfilename('fullpath'),'.m'];
%     saveSequenceCopy(fpathfull,sq.directory,varargin);
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
