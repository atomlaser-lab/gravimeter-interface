function sq = makeBEC(varargin)

opt = SequenceOptions;

if nargin == 1
    %
    % If first argument is of type GravimeterOptions, use that
    %
    if ~isa(varargin{1},'SequenceOptions')
        error('If using one argument it must be of type SequenceOptions');
    end
    opt = opt.replace(varargin{1}); 
elseif mod(nargin,2) == 0
    opt.set(varargin{:});
elseif mod(nargin - 1,2) == 0 && isa(varargin{1},'SequenceOptions')
    opt = opt.replace(varargin{1});
    opt.set(varargin{2:end});    
else 
    error('You must supply either a SequenceOptions object, a set of name/value pairs for options, or a SequenceOptions object followed by name/value pairs');
end

%% Create a conversion object to handle conversions to volts
convert = RunConversions;

%% Initialize sequence - defaults should be handled here
sq = initSequence;

%    U/D bias field converter to amps to volts (possible values are -0.58 A to 14 A) which corresponds to a voltage from (2.823V to -0.14V)
UD = @(x) x* -0.2031 + 2.707; %this converts the input value in amps to the appropriate voltage
% UDreverse = @(y) y*-4.924 + 13.33; %This is just the inverted function to get the amps out of the voltage in case needed.

sq.find('Imaging Freq').set(convert.imaging(opt.detuning));
sq.find('3D MOT Freq').set(convert.mot_freq(-23));    %Use -25 MHz for 4 s loading times, -27.5 MHz for 6 s
sq.find('Repump freq').set(convert.repump_freq(-2));
sq.find('50w ttl').set(1);
sq.find('25w ttl').set(1);
sq.find('50w amp').set(convert.dipole50(15));
sq.find('25w amp').set(convert.dipole25(15)); 

%% Set up the MOT loading values                
sq.find('MOT coil TTL').set(1);     %Turn on the MOT coils
sq.find('3d coils').set(convert.mot_coil(2.58));
sq.find('bias u/d').set(0);
sq.find('bias e/w').set(0);
sq.find('bias n/s').set(0);

Tmot = 4;                           %6 s MOT loading time
sq.delay(Tmot);                     %Wait for Tmot
%% Compressed MOT stage
%Turn off the 2D MOT and push beam 10 ms before the CMOT stage
sq.find('2D MOT Amp TTL').before(10e-3,0);
sq.find('push amp ttl').before(10e-3,0);

%Increase the cooling and repump detunings to reduce re-radiation
%pressure, and weaken the trap
sq.find('3D MOT freq').set(convert.mot_freq(-25));
sq.find('repump freq').set(convert.repump_freq(-6));
sq.find('3D coils').set(convert.mot_coil(1.38));
sq.find('bias e/w').set(5*0);
sq.find('bias n/s').set(7*0);
sq.find('bias u/d').set(1.75*0);

Tcmot = 12.5e-3;                      %10 ms CMOT stage
sq.delay(Tcmot);                    %Wait for time Tcmot
%% PGC stage
Tpgc = 25e-3;
%Define a function giving a 100 point smoothly varying curve
t = linspace(0,Tpgc,50);
f = @(vi,vf) sq.linramp(t,vi,vf);

%Smooth ramps for these parameters
sq.find('3D MOT Amp').after(t,f(5,3.8));
sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-60)));
sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),convert.mot_coil(0.45)));
sq.find('repump freq').set(convert.repump_freq(-9.25));

sq.delay(Tpgc);
% Turn off the repump field for optical pumping - 1 ms
T = 1e-3;
sq.find('repump amp ttl').set(0);
% sq.find('Top repump shutter').set(1);
sq.find('liquid crystal repump').set(-2.22);
sq.find('bias u/d').set(0);
sq.find('bias e/w').set(0);
sq.find('bias n/s').set(7.5*0);
sq.delay(T);

% Load into magnetic trap
sq.find('liquid crystal bragg').set(-3);
sq.find('3D mot amp ttl').set(0);
sq.find('MOT coil ttl').set(1);
sq.find('3D coils').set(convert.mot_coil(6));
% sq.delay(2.50);

% %% Microwave evaporation
% %
% % Provide detunings in MHz from the Rb hyperfine splitting
% %
sq.delay(20e-3);
evapRate = 15;
evapStart = 55;
evapEnd = 10;
Tevap = (evapStart-evapEnd)/evapRate;
t = linspace(0,Tevap,100);
sq.find('mw freq').after(t,convert.microwave(sq.linramp(t,evapStart,evapEnd)));
sq.find('mw amp ttl').set(1);
sq.delay(Tevap);
sq.find('mw amp ttl').set(0);
% sq.delay(3);
%% Weaken trap while MW frequency fixed
% Trampcoils = 180e-3;
% t = linspace(0,Trampcoils,50);
% sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),convert.mot_coil(1)));
% sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),0));
% sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),0));
% sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),-0.12*0));
% sq.delay(Trampcoils);

% Trampcoils = 0.5/2;
% t = linspace(0,Trampcoils,50);
% sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
% sq.find('bias e/w').after(t,sq.linramp(t,sq.find('bias e/w').values(end),0));
% sq.find('bias n/s').after(t,sq.linramp(t,sq.find('bias n/s').values(end),0));
% sq.find('bias u/d').after(t,sq.linramp(t,sq.find('bias u/d').values(end),0));
% sq.delay(Trampcoils);
% sq.find('mot coil ttl').set(0);

% % Optical evaporation
% 
% % Ramp down magnetic trap in 1 s
% 
Trampcoils = 1;
t = linspace(0,Trampcoils,101);
sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
sq.find('mw amp ttl').anchor(sq.find('3d coils').last).before(100e-3,0);
sq.find('mot coil ttl').at(sq.find('3d coils').last,0);
% 
% % At the same time, start optical evaporation
% 
sq.delay(30e-3);
Tevap = 2;
t = linspace(0,Tevap,200);
sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(opt.dipoles),0.4));
sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(opt.dipoles),0.4));

sq.find('bias e/w').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias e/w').values(end),10));
sq.find('bias n/s').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias n/s').values(end),0));
sq.find('bias u/d').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias u/d').values(end),0));
sq.delay(Tevap);

% sq.delay(opt.params(1));



end