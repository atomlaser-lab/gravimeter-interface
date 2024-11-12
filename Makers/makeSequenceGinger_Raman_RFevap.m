function varargout = makeSequenceGinger(varargin)
%% Parse input arguments
opt = SequenceOptions('load_time',15,'detuning',0,'tof',20e-3,'redpower',2,...
    'raycus',2);

% check if SequenceOptions already exists in the workspace
if ~exist('SequenceOptions','class')
    % create a new instance of SequenceOptions if it does not exist
    opt = SequenceOptions('load_time',5,'detuning',0,'tof',17.3e-3,'redpower',2,'raycus',2);
    %     assignin('caller', 'opt', opt);
else
    % use the existing instance of SequenceOptions
    opt = SequenceOptions();
end

% check input arguments and update options accordingly
if nargin == 0
    if evalin('base', 'exist(''opt'', ''var'')') == 1
        % use the existing instance of SequenceOptions
        opt = evalin('base', 'opt');
    else
        opt = SequenceOptions('load_time', 5, 'detuning', 0, 'tof', 17.3e-3, 'redpower', 2, 'raycus', 2);
        assignin('base', 'opt', opt);
    end

elseif nargin == 1
    if ~isa(varargin{1},'SequenceOptions')
        error('If using only one argument it must of type SequenceOptions');
    end
    opt.replace(varargin{1});
elseif mod(nargin,2) == 0
    opt.set(varargin{:});
elseif mod(nargin - 1,2) == 0 && isa(varargin{1},'SequenceOptions')
    opt.replace(varargin{1});
    opt.set(varargin{2:end});
    assignin('base', 'opt', opt);

else
    error('Either supply a single SequenceOptions argument, or supply a set of name/value pairs, or supply a SequenceOptions argument followed by name/value pairs');
end

if nargout == 0
    r.make(opt).urun(@Abs_Analysis_Fancy);
    assignin('base', 'opt', opt);
else
    varargout{1} = opt;
end

% opt = ensureSequenceOptionsAvailable();

ImageFreq = opt.detuning*0.6238/6 + 8.3; %assuming imaging field is set to 5V
dipole_field = 5;
imaging_field = dipole_field; %dipole_field
ImageAmp = 8; %5 until and including mag load, 9 for high OD samples on horizontal imaging % 5 for thermal

%make SLM pattern
l = 0; %LG charge
phase = 90;
% phase=opt.params;
sign = [-1 -1 -1];
f = 0; %79 for charge 1 %focal length
number_of_pulses = 1;
dir = 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\Image Files\512\test\'; %to save images to be sent to the SLM
LG_grating_generator(l,phase,sign,f,dir,number_of_pulses);

%% Set Sequence

CMOT = 1;
PGC = 1;
Mag_Trap = 1;
RF_Knife = 1;
ODT_Load = 1;
VMG_Setup = 0;
ODT_Evap = 0;
ARP = 0;
RF_Pi_Pulse = 0; %0
VMG_Run = 1;
SG_Pulse = 0; %0

%% Initialize sequence
sq = initSequence;  %load default values (OLD MOT values are default)
sq.find('87 imag freq').set(8.3);
sq.find('87 imag amp').set(8);

%% MOT loading
%
% We use CD channel 0b10 = 2 for loading the MOT
%

% Tmot = opt.load_time;
Tmot = 3; %15 for vertical

% sq.find('C6 - N/C').set(1).after(100e-6,0); %SLM trigger

sq.find('2DMOT').set(1);
sq.find('3DMOT').set(1);
sq.find('87 push').set(1);
% 3D MOT beam settings
sq.find('3DMOT Freq').set(FtoV('trap',26));
sq.find('3DMOT amp').set(TrapPtoV('trap',0.7)); % WAS: 0.4 NOW 0.7
% 3D repump beam settings
sq.find('87 repump').set(1);
sq.find('Repump Switch').set(0);
sq.find('87 repump freq').set(FtoV('repump',0));
sq.find('87 repump amp').set(TrapPtoV('repump',1));
% 3D coil settings
sq.find('H-Bridge Quad').set(1);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(1);
sq.find('CD2').set(dBtoV('normal',24)); %Coarse control of 3D coils
sq.find('CD Fine/Fast').set(dBtoV('fine',0)); % fine control of 3D coils
% Bias coil settings
sq.find('Earth Bias 1').set(3); %6
sq.find('Earth Bias 2').set(5); %4
% sq.find('Earth Bias 3').set(0); %2

% sq.find('LG Shutter').set(1).after(1,0); %triggering shutter

%Delay for the load time
sq.delay(Tmot);

%
% Turn off the 2D MOT and coils as well as the push beam
%
sq.find('2D MOT Coils').before(10e-3,0);
sq.find('2DMOT').before(10e-3,0);
sq.find('87 push').before(10e-3,0);
sq.find('85 push').before(10e-3,0);

sq.find('C6 - N/C').set(1).after(100e-6,0); %trigger SLM

%%% CMOT sequence %%%
if(CMOT)
    %
    % Apply a compressed MOT sequence to temporarily increase the density by
    % reducing spontaneous emission.  We switch to CD channel 0b00 = 0 because
    % it is the fast channel
    %

    Tcmot = 30e-3;
    t = 0:2e-3:Tcmot; %%NEW

    %3D Coils
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    sq.find('CD0 Fast').set(0);
    sq.find('CD0 Fast').set(dBtoV('normal',0)); %%NEW
    sq.find('CD Fine/Fast').set(dBtoV('fine',10));

    %Trapping light
    sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),FtoV('trap',53.4))); % WAS 53.4
    sq.find('3DMOT amp').set(TrapPtoV('trap',1.45)); % WAS 1.0 NOW 1.45

    %Repump
    sq.find('87 repump freq').set(FtoV('repump',12)); %-7
    sq.find('87 repump amp').set(TrapPtoV('repump',0.05)); %0.05

    sq.delay(Tcmot);

    %%% PGC sequence %%%
    if(PGC)
        %
        % Apply polarization gradient cooling to reduce the temperature of the
        % atoms.  We use CD channel 0b00 = 0 as it is the fast-switching channel
        %

        Tpgc = 6e-3;
        t = linspace(0,Tpgc,10); %%NEW

        %Earth Coils
        sq.find('Earth Bias 1').set(1.2); %Probably N/S
        sq.find('Earth Bias 2').set(6.5);%Pr0bably U/D
        %         sq.find('Earth Bias 3').set(0.05); % Probably E/W

        %3D Coils
        sq.find('CD fine/fast').set(0);
        sq.find('CD0 Fast').set(0);

        %Trapping light
        sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),FtoV('trap',79.8)));
        sq.find('3DMOT amp').after(t,sq.minjerk(t,sq.find('3DMOT amp').values(end),TrapPtoV('trap',0.85))); % WAS 0.8, NOW 0.85

        %Repump
        sq.find('87 repump freq').set(FtoV('repump',11));%-4.8
        sq.find('87 repump amp').set(TrapPtoV('repump',0.05));%0.004

        sq.delay(Tpgc);

        %
        % Turn off repump field so that atoms are optically pumped into the F = 1
        % manifold.
        %

        Tdepump = 1e-3;
        sq.find('repump switch').set(1); %fiber switch off (it's inverted)
        sq.find('87 repump').set(0);
        sq.find('87 repump amp').set(0);
        sq.find('85 repump').set(0);
        sq.delay(Tdepump);
        sq.find('3DMOT').set(0);

        %%% Load into magnetic trap %%%
        if(Mag_Trap)
            % Load into the magnetic trap at a high gradient.  We switch quickly to a
            % low value and then ramp up to the target value
            %

            sq.find('Earth Bias 1').set(0);
            sq.find('Earth Bias 2').set(0);
            %             sq.find('Earth Bias 3').set(0);

            Tmagload = 150e-3;
            t = 0:10e-3:Tmagload;
            dBLoad = 110;
            sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',dBLoad/2),dBtoV('normal',dBLoad)));
            sq.find('CD Fine/Fast').set(dBtoV('fine',0));
            %             sq.delay(Tmagload);

            Toptload = 400e-3;
            t = 0:20e-3:Toptload;
            sq.find('Raycus TTL').set(1);
            sq.find('Redpower TTL').set(1);
            sq.find('Raycus CW').after(t,sq.minjerk(t,0,DipolePtoV('raycus',8))); %8
            sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower',13))); %13
            sq.find('MOT bias').set(1);
            sq.find('MOT bias coil').after(t,sq.linramp(t,0,dipole_field));
            sq.delay(max(Toptload,Tmagload));

            % sq.delay(1);

            %%% RF evaporation %%%
            if (RF_Knife)
                %
                % Remove hot atoms from the sample using RF transitions between the trapped
                % |F = 1, m_F = -1> state and the untrapped |F = 1, m_F = 0> state.  All
                % frequencies are in MHz
                %
                rf_start = 16;
                rf_end = 10; %1
                rf_rate = 3; %MHz/s %3.5
                Tevap = (rf_start - rf_end)/rf_rate;
                rf_ramp_type = 'lin';
                rf_exp_time_constant = 2;
                t = linspace(0,Tevap,50);

                sq.find('RF atten').set(1);
                if strcmpi(rf_ramp_type,'exp')
                    sq.find('RF frequency').after(t,sq.expramp(t,RFtoV(rf_start),RFtoV(rf_end),rf_exp_time_constant)); %ramp rf frequency from 4 to -2.667
                elseif strcmpi(rf_ramp_type,'lin')
                    sq.find('RF frequency').after(t,sq.linramp(t,RFtoV(rf_start),RFtoV(rf_end)));
                end
                sq.delay(Tevap);

                sq.find('RF atten').set(0);
                sq.find('RF Frequency').set(RFtoV(20));
                %
                %                                                 % Turn off magnetic trap (use for dipole alignment)
                %                                                 sq.find('CD0 Fast').set(0);
                %                                                 sq.find('CD Fine/Fast').set(0);
                %                                                 sq.delay(10e-3 - opt.tof);

                %%% Ramp down coils %%%
                if(ODT_Load)

                    %for main coils
                    sq.find('CD0 Fast').set(0);
                    sq.find('CD Fine/Fast').set(0);

                    %                                                             %                     to hold in dipole trap
                    %                                                             sq.delay(1);

                    %%% VMG Setup %%%
                    if (VMG_Setup)
                        %
                        %
                        % Switch to Helmholtz configuration for state preparation and blow away F=2
                        %
                        %

                        sq.find('CD bit 0').set(0);
                        sq.find('CD bit 1').set(0);

                        sq.find('H-Bridge Quad').set(0);
                        sq.delay(50e-6);
                        sq.find('H-Bridge Helm').set(1);
                        sq.delay(100e-6);
                        Tramp = 1; %1
                        t = linspace(0,Tramp,51);
                        sq.find('CD0 Fast').after(t,sq.minjerk(t,0,dBtoV('normal',20))); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCOMMENT ME FOR PROPER STATE PREP
                        %sq.find('MOT bias coil').after(t,sq.minjerk(t,sq.find('MOT bias coil').values(end),0));
                        sq.delay(Tramp);

                    end

                    %%% Optical evaporation %%%
                    if (ODT_Evap)
                        Tevap = 3.5;
                        TC = 0.5; %0.5
                        t = linspace(0,Tevap,150);
                        final_dipole.RP = 2; %2 %2.5 %1.5 %1.3 | WAS 1.0
                        final_dipole.FA = 2; %2 %2.5 %1.5 %0.9 | WAS 1.0
                        %                                                                         final_dipole.RP = opt.raycus;
                        %                                                                         final_dipole.FA = opt.raycus;
                        sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP),TC));
                        sq.find('Raycus CW').after(t,sq.expramp(t,sq.find('Raycus CW').values(end),DipolePtoV('raycus',final_dipole.FA),TC));
                        sq.delay(Tevap);
                    end

                                        %blow away F=2 atoms
                                        Tblow = 10e-6;
                                        sq.find('3DMOT').set(1);
                                        sq.find('variable wave plate').set(-2.42);
                                        sq.find('3DMOT Freq').set(7.65); %FtoV('trap',0)
                                        sq.find('3DMOT amp').set(TrapPtoV('trap',1)); % WAS: 0.4 NOW 0.7
                                        sq.delay(Tblow);
                                        sq.find('3DMOT').set(0);
                    
                                        %
                                        % Trigger the DDS
                                        %
                                        sq.ddsTrigDelay = sq.time;
                                        sq.find('DDS TTL').before(10e-3,1).after(10e-3,0);%.after(1e-3,1);

                    %%% ARP with DDS %%%
                    if (ARP)

%                         %blow away F=2 atoms
%                         Tblow = 10e-6;
%                         sq.find('3DMOT').set(1);
%                         sq.find('variable wave plate').set(-2.42);
%                         sq.find('3DMOT Freq').set(7.65); %FtoV('trap',0)
%                         sq.find('3DMOT amp').set(TrapPtoV('trap',1)); % WAS: 0.4 NOW 0.7
%                         sq.delay(Tblow);
%                         sq.find('3DMOT').set(0);
% 
%                         %
%                         % Trigger the DDS
%                         %
%                         sq.ddsTrigDelay = sq.time;
%                         sq.find('DDS TTL').before(10e-3,1).after(10e-3,0);%.after(1e-3,1);

                        sq.dds(1).set(38,0,0);
                        sq.dds(2).set(110,0,0);
                        sq.delay(10e-6);
                        sq.find('RF Switch').set(1);
                        Tarp = 10e-3; %25e-3 %80e-3
                        t = linspace(0,Tarp,501);
                        %df = 38.35 31/10/22
                        df = 38.5 + 3*sq.linramp(t,-0.5,0.5);
                        w = Tarp/3;
                        % amp = 0.075*sech((t - Tarp/2)/w).^2;
                        % amp = 1.3e-3*sech((t - Tarp/2)/w).^2;
                        %                             amp = 1.5e-3*sech((t - Tarp/2)/w).^2;
                        %                             amp = 0.05*0.022*sech((t - Tarp/2)/w).^2;
                        %                             amp = 0.001*ones(size(t));
                        amp = 5e-3*exp(-2*(t - Tarp/2).^2/w^2);
                        sq.dds(1).after(t,df,amp,0);
                        sq.dds(2).after(t,110,0,0);
                        sq.delay(Tarp);
                        sq.dds(1).set(110,0,0);
                        sq.dds(2).set(110,0,0);
                        sq.find('RF Switch').set(0);
                    end

                    %%% RF Pi Pulse from |1,-1> to |1,0>
                    if(RF_Pi_Pulse)
                        pulse_freq = 38; %38.263 was quoted to work on 16/9/22
                        sq.dds(1).set(pulse_freq,0,0);
                        sq.dds(2).set(110,0,0);
                        sq.delay(50e-6);
                        sq.find('RF Switch').set(1);
                        sq.dds(1).set(pulse_freq,0.075*1,0);
                        sq.dds(2).set(110,0,0);
                        sq.delay(10e-6);
                        sq.dds(1).set(pulse_freq,0,0);
                        sq.dds(2).set(110,0,0);
                        sq.find('RF Switch').set(0);
                    end


                end
            end
        end

    end
end
%% Drop atoms
% sq.delay(1);
% sq.find('Earth Bias 1').set(0);
% sq.find('Earth Bias 2').set(0);
% sq.find('Earth Bias 3').set(0);

timeAtDrop = sq.time;
% sq.find('Probe').set(1).after(1e-3,0);
%
% This trigger delay is necessary because the DDS instructions start when
% the DDS trigger occurs
%
% sq.ddsTrigDelay = timeAtDrop;
% sq.find('DDS TTL').before(10e-3,1).after(10e-3,0).after(1e-3,1);
%
% Set all other channels to 0

sq.find('2D MOT Coils').set(0);
sq.find('3DMOT').set(0);
sq.find('87 repump').set(0);
sq.find('87 repump amp').set(0);
sq.find('CD0 Fast').set(dBtoV('normal',0)); %REMEMBER TO TURN BACK TO 2 WHEN USING RAMAN??
sq.find('CD2').set(0);
sq.find('CD Fine/Fast').set(0);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
sq.find('Redpower CW').set(0);
sq.find('Redpower TTL').after(100e-6,0);
sq.find('Raycus CW').set(0);
sq.find('Raycus TTL').after(5e-3,0);
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));

%marker for bottom of cell at 45ms from drop (22ms for horizontal imaging)
% sq.find('C6 - N/C').set(0).after(45e-3,1).after(10e-6,0);

% % sq.delay(100e-6);
% % sq.find('H-bridge helm').set(0);
% % sq.delay(50e-6);
% % sq.find('H-bridge quad').set(1);

%%% VMG! %%%
if(VMG_Run)
    sq.anchor(timeAtDrop);
    raman_delay_drop = 5e-3; %5e-3 %10.8e-3
    interrogation_time = 1e-3; %6e-3
    Traman_pi = opt.params; %18e-6
    F2_tof = 7e-3; %20.1e-3 %11e-3
    F1_tof = F2_tof + 6e-3; %3e-3 minimum with smallest frame for BEC
    F1_minus_F2_tof = F1_tof - F2_tof; %needs to be 12ms minimum for full frame
    AI_phi_G = 0; %abs(180 - opt.params*180*1e3)
    AI_phi_LG = 0;
    sq.anchor(timeAtDrop + raman_delay_drop);
    Power_raman_G = 1;  %0.3 %sideband
    Power_raman_LG = 0.1; %1 %carrier
    Delta_raman = 19.55; %20.02 + 1.790493109783823e-2 + 910.6846/1e6

    % state prep

%     % %     %small field with the MOT coils
%     sq.find('H-Bridge Quad').set(0);
%     sq.delay(50e-6);
%     sq.find('H-Bridge Helm').set(1);
%     sq.delay(100e-6);
%     sq.find('CD0 Fast').before(raman_delay_drop,0.2); %7.52G
    sq.find('MOT bias').before(raman_delay_drop+1,0);
%     %      sq.find('LG Shutter').before(raman_delay_drop+1,1).after(F2_tof+raman_delay_drop+1,0); %opening shutter
%     % %                 sq.delay(1e-3);

    %     single Raman pulse
    sq.dds(1).set(110+Delta_raman/4,Power_raman_G,0);
    sq.dds(2).set(110-Delta_raman/4,Power_raman_LG,0);
    sq.delay(Traman_pi);
    sq.dds(1).set(110,0,0);
    sq.dds(2).set(110,0,0);
    sq.delay(10e-6);

    %         %     BS1
    %         sq.dds(1).set(110+Delta_raman/4,Power_raman_G,0); %sideband
    %         sq.dds(2).set(110-Delta_raman/4,Power_raman_LG,0); %carrier
    %         sq.delay(Traman_pi/2);
    %         sq.dds(1).set(110+Delta_raman/4,0,AI_phi_G);
    %         sq.dds(2).set(110-Delta_raman/4,0,AI_phi_LG);
    %     %     sq.find('LG Shutter').set(1).after(100e-6,0); %trigger SLM
    %
    %         sq.delay(interrogation_time);
    %
    %         %     BS2
    %         sq.dds(1).set(110+Delta_raman/4,Power_raman_G,AI_phi_G); %sideband
    %         sq.dds(2).set(110-Delta_raman/4,Power_raman_LG,AI_phi_LG); %carrier
    %         sq.delay(Traman_pi/2);
    %         sq.dds(1).set(110+Delta_raman/4,0,0);
    %         sq.dds(2).set(110-Delta_raman/4,0,0);
    %         sq.delay(10e-6);

    %turning field off and imaging field on
    sq.find('CD0 Fast').set(0); %0.133 ~0.5G
    sq.find('MOT bias').set(1);
    %             sq.delay(raman_delay_drop + Traman_pi);

    % blow away atoms in F=2
    % sq.delay(0.5e-3);
    % sq.find('87 imag').set(1);
    % sq.delay(1e-3);
    % sq.find('87 imag').set(0);

end

%%% S -G Field (SG pulse) %%%
% SG_Pulse = 0;
if(SG_Pulse)

    % HEY RYAN! MAKE SURE THIS IS UNCOMMENTED TO SWITCH BACK TO ANTI-HELMHOLTZ CONFIGURATION FOR STERN-GERLACH
    sq.delay(10e-6);
    sq.find('CD0 Fast').set(0);
    sq.delay(100e-6);
    sq.find('H-bridge helm').set(0);
    sq.delay(50e-6);
    sq.find('H-bridge quad').set(1);

    sq.anchor(timeAtDrop);
    sq.delay(raman_delay_drop+Traman_pi+interrogation_time+0.3e-3); %1e-3
    SG_Pulse_time = 3.3e-3; %3.45 %2.45e-3 %2.9e-3 %8e-3
    SG_Amp = 70; %25
    t = linspace(0,SG_Pulse_time,50);
    sq.find('CD0 Fast').after(t,sq.minjerk(t,0,dBtoV('normal',SG_Amp)));
    sq.delay(SG_Pulse_time);
    sq.find('CD0 Fast').after(t,sq.minjerk(t,sq.find('CD0 Fast').values(end),0));

    % F2_tof = 3e-3;
    % F1_tof = 15e-3;
end
%% Take Absorption Image

sq.anchor(timeAtDrop);
sq.camDelay = timeAtDrop - 2;

% makeImagingSequence(sq,'tof',opt.tof,'pulse time',80e-6,'repump delay',100e-6,...
%     'repump time',200e-6,'cam time',5e-6,'cycle time',100e-3,...
%     'manifold',2,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
%     'fibre switch delay',1e-3,'imaging_field',imaging_field,'image type','horizontal');

makeImagingSequence_4_Images(sq,'tof',F2_tof,'tof2',F1_minus_F2_tof,'pulse time',80e-6,'repump delay',100e-6,...
    'repump time',200e-6,'cam time',50e-6,'cycle time',100e-3,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
    'fibre switch delay',1e-3,'imaging_field',imaging_field,'image type','horizontal');

% turn off the dipoles
sq.find('Redpower CW').set(0);
sq.find('Redpower TTL').after(100e-6,0);
sq.find('Raycus CW').set(0);
sq.find('Raycus TTL').after(5e-3,0);
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));
sq.find('3DMOT').set(0);
sq.find('87 repump amp').set(0);
sq.find('CD0 Fast').set(0);
% sq.find('C6 - N/C').set(1).after(100e-6,0).after(10e-3,1).after(100e-6,0);
% sq.delay(100e-6);
% sq.find('H-bridge helm').set(0);
% sq.delay(50e-6);
% sq.find('H-bridge quad').set(1);

sq.waitFromLatest(0.25);
setSafeValues(sq);
% sq.delay(5);

if nargout == 0
    r = RemoteControl;
    r.upload(sq.compile);
    r.run;
else
    varargout{1} = sq;
end

end