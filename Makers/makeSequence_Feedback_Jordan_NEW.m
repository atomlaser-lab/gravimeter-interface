function varargout = makeSequence_Feedback_Jordan_NEW(varargin)
%% BEC_BUILDER *Make BEC great again!*
%% 
% 
%% 
% 
% Check if any variable is an instance of SequenceOptions
% Get all variable names in the workspace
allVarNames = evalin('base', 'who');
% Check if any variable is an instance of SequenceOptions
seqOptVarIndices = cellfun(@(varName) isa(evalin('base', varName), 'SequenceOptions'), allVarNames);
seqOptVarNames = allVarNames(seqOptVarIndices);
if isempty(seqOptVarNames)
    if  nargin ==1 && ~exist(varargin{1},'var')
        error('define the sequence option variable name first');
    end
    % Create a new instance of SequenceOptions if it does not exist
    opt = SequenceOptions('load_time', 5, 'detuning', 0, 'tof', 17.3e-3, 'redpower', 2, 'keopsys', 2);
    cprintf('Keywords','Initialising Sequence Options\n')
    assignin('base', 'opt', opt);
else
    if numel(seqOptVarNames) > 1
        error('Multiple SequenceOptions found in the workspace');
    end
    % Use the existing instance of SequenceOptions
    opt = evalin('base', seqOptVarNames{1});
end
% Check input arguments and update options accordingly
if nargin == 1
    if ~isa(varargin{1}, 'SequenceOptions')
        error('If using only one argument, it must be of type SequenceOptions');
    end
elseif nargin > 1
    if ~isa(varargin{1}, 'SequenceOptions')
        error('First argument must be of type SequenceOptions');
    end
    opt.replace(varargin{1});
    opt.set(varargin{2:end});
end
if nargout == 0
    % Search for the opt variable in the base workspace
    varExists = evalin('base', 'exist(''opt'', ''var'')');
    if varExists
        opt = evalin('base', 'opt');
    end
else
    varargout{1} = opt;
end
%% 
% 
%% Camera Type
% Use to select the camera in-use
Camera = "low res";
% Drop_time = evalin('base', 'opt.tof');
%% 
% 
%% Cooling sequence
% Not if stage _j_ is off, all stages after _j_ are off. 
CMOT = 1  ;
PGC = 1  ;
Load_Mag_Evap = 1  ;
RF_Evaporation = 1  ;
Load_Opt_Evap = 1  ;
Optical_Evaporation = 1  ;
% Final_dipole-common_power = 2;
%% 
%% Display Information Status
Display_info = 1 ;
%% 
%% Initialisation of useful conversion for Dipoles, Imaging...
% ImageFreq = opt.detuning*0.6238/6 + 8.3; %assuming imaging field is set to 5V
% ImageFreq = opt.params; %assuming imaging field is set to 5V
ImageFreq = (opt.detuning+43.5)^4*6e-8-(opt.detuning+43.5)^3*8e-7-(opt.detuning+43.5)^2*0.0004+(opt.detuning+43.5)*0.0796+5.4703; %assuming imaging field is set to 5V
% ImageFreq = opt.detuning^4*6e-8-opt.detuning^3*8e-7-opt.detuning^2*0.0004+opt.detuning*0.0796+5.7403; %assuming imaging field is set to 5V
dipole_field = 5; %
imaging_field = dipole_field; %5
ImageAmp = 9; %5 until and including mag load, 7 for high OD samples on horizontal imaging
%% 
%% Initialistaiont of the MOT sequence
% sq = initSequence;  %load default values
% sq.find('87 imag freq').set(8.3); % 8.35
% sq.find('87 imag amp').set(8); % 8
%
% sq.find('2DMOT').set(1);
% sq.find('3DMOT').set(1);
% sq.find('87 push').set(1);
%
% % 3D MOT beam settings
% % sq.find('3DMOT Freq').set(FtoV('trap',16));
% sq.find('3DMOT Freq').set(7.25); % 6.7
% % sq.find('3DMOT amp').set(TrapPtoV('trap',1));
% sq.find('3DMOT amp').set(7); %8
%
% % 3D repump beam settings
% sq.find('Repump Switch').set(0);
% sq.find('87 repump').set(1);
% % sq.find('87 repump freq').set(FtoV('repump',3)); %0
% sq.find('87 repump freq').set(4.68); %4.68
%
% % sq.find('87 repump amp').set(TrapPtoV('repump',0.93)); %1
% sq.find('87 repump amp').set(8);
%
% % 3D coil settings
% sq.find('H-Bridge Quad').set(1);
% sq.find('CD bit 0').set(0);
% sq.find('CD bit 1').set(1);
% % sq.find('CD2').set(dBtoV('normal',14)); %Coarse control of 3D coils %11
% sq.find('CD2').set(1.18); %1.8
% % sq.find('CD Fine/Fast').set(dBtoV('fine',8)); %8 fine control of 3D coils
% sq.find('CD Fine/Fast').set(0); %8 fine control of 3D coils
%
% % sq.find('85 repump amp').set(5.5); % 3DMOT bias
% % sq.find('MOT bias').set(1); %switching imaging coils on
% % sq.find('MOT bias coil').set(1.5); %imaging coils
%
% % sq.find('Earth Bias 1').set(10);
% % sq.find('Earth Bias 2').set(2);
% % sq.find('Earth Bias 3').set(8);
%
% % Tmot = 7;
% Tmot = opt.load_time;
% sq.delay(Tmot);
sq = initSequence; %load default values
sq.find('87 imag freq').set(8.3); % 8.35
sq.find('87 imag amp').set(8); % 8
sq.find('2DMOT').set(1);
sq.find('3DMOT').set(1);
sq.find('87 push').set(1);
% 3D MOT beam settings
% sq.find('3DMOT Freq').set(FtoV('trap',16));
%% sq.find('3DMOT Freq').set(7.250); % 6.7 %*was 6.6*
sq.find('3DMOT Freq').set(6.6); % 6.7 %*was 6.6*
% sq.find('3DMOT amp').set(TrapPtoV('trap',1));
%% sq.find('3DMOT amp').set(7.000); %8 %*was 4.4*
sq.find('3DMOT amp').set(4.4); %8 %*was 4.4*
% 3D repump beam settings
sq.find('Repump Switch').set(0);
sq.find('87 repump').set(1);
% sq.find('87 repump freq').set(FtoV('repump',3)); %0
%% sq.find('87 repump freq').set(4.680); %4.68 %*was 4.4*
sq.find('87 repump freq').set(4.4); %4.68 %*was 4.4*
% sq.find('87 repump amp').set(TrapPtoV('repump',0.93)); %1
%% sq.find('87 repump amp').set(8); %*was 7*
sq.find('87 repump amp').set(7); %*was 7*
% 3D coil settings
sq.find('H-Bridge Quad').set(1);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(1);
% sq.find('CD2').set(dBtoV('normal',14)); %Coarse control of 3D coils %11
%% sq.find('CD2').set(1.180); %1.8 %*was 1.8*
sq.find('CD2').set(1.80); %1.8 %*was 1.8*
% sq.find('CD Fine/Fast').set(dBtoV('fine',8)); %8 fine control of 3D coils
sq.find('CD Fine/Fast').set(0); %8 fine control of 3D coils
% sq.find('85 repump amp').set(5.5); % 3DMOT bias
% sq.find('MOT bias').set(1); %switching imaging coils on
% sq.find('MOT bias coil').set(1.5); %imaging coils
sq.find('Earth Bias 1').set(3); %6
sq.find('Earth Bias 2').set(5); %4
sq.find('Earth Bias 3').set(0); %2
% Tmot = opt.params;
Tmot = opt.load_time;
sq.delay(Tmot);
%% 
% 
%% Compression MOT sequence
if CMOT == 1
%% 
% 
% 
% *Turn off the 2D MOT and coils as well as the push beam*
    sq.find('2D MOT Coils').before(10e-3,0);
    sq.find('2DMOT').before(10e-3,0);
    sq.find('87 push').before(10e-3,0);
    sq.find('85 push').before(10e-3,0);
%% 
% 
    %     Tcmot = 0.012; %15e-3
    Tcmot = 30e-3; %10e-3
%% 
% 
    t = 0:1e-3:Tcmot;
    %     %3D Coils (old)
    %     sq.find('CD bit 0').set(0);
    %     sq.find('CD bit 1').set(1);
    %     sq.find('CD2').set(.6); %.6
    %     sq.find('CD Fine/Fast').set(0);
    %3D Coils (new)
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    sq.find('CD0 Fast').set(dBtoV('normal',0));
    %     sq.find('CD Fine/Fast').set(dBtoV('fine',5));
    sq.find('CD Fine/Fast').set(7.5); %2
    %Trapping light
    %     sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),FtoV('trap',35))); %45
    sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),4.6)); %4.5
    sq.find('3DMOT amp').set(6); %4.5
    %Repump
    %     sq.find('87 repump freq').set(FtoV('repump',6)); %-7.5 %-4
    sq.find('87 repump freq').set(6.25); %6.2
    sq.find('87 repump amp').set(6); %5
    sq.delay(Tcmot);
%% 
% *Back to contents*
%% PGC
    if PGC == 1
%% 
% 
        %         Tpgc = 0.005; %usually at 5 ms...
        Tpgc = 6e-3; %4e-3
%% 
% 
        t = 0:1e-3:Tpgc;
        % %
        sq.find('Earth Bias 1').set(1.2);  %Probably N/S
        sq.find('Earth Bias 2').set(6.5);  %Pr0bably U/D
        sq.find('Earth Bias 3').set(0.05); % Probably E/W
%         sq.find('CD fine/fast').set(0.2); %0.2
        sq.find('CD fine/fast').set(0); %0.2
        sq.find('CD0 Fast').set(0); %1
%         sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),FtoV('trap',77.5)));
        sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),2.3));
        sq.find('3DMOT amp').after(t,sq.minjerk(t,sq.find('3DMOT amp').values(end),4.1)); %4
        %         sq.find('87 repump freq').set(FtoV('repump',11));
        sq.find('87 repump freq').set(2.3); %2.25
        sq.find('87 repump amp').set(3); %3
        sq.delay(Tpgc);
% Optical pumping into F =1 manifold
        Tdepump = 1e-3;
        sq.find('repump switch').set(1); %fiber switch off (it's inverted)
        sq.find('87 repump').set(0);
        sq.find('87 repump amp').set(0);
        sq.find('85 repump').set(0);
        sq.delay(Tdepump);
        sq.find('3DMOT').set(0);
%% 
% *Back to contents*
%% 
%% Load into the Magnetic Trap
        if Load_Mag_Evap == 1
            %             Tmagload = 0; %150e-3
            Tmagload = 200e-3; %140e-3
%% 
% 
                       sq.find('Earth Bias 1').set(0);
                       sq.find('Earth Bias 2').set(0);
                       sq.find('Earth Bias 3').set(0);
            % %
            %             t = linspace(0,Tmagload,50);
            %             dBLoad = 20; %110 Ryan has that
            %             sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',dBLoad/2),dBtoV('normal',dBLoad)));
            %             %             sq.find('CD Fine/Fast').set(dBtoV('fine',0)); %Ryan Has 0
            %             sq.find('CD Fine/Fast').set(0);
            %
            %
            %             sq.delay(Tmagload); %Tmagload
            %
            %             Toptload = 400e-3; %400e-3
            %             t = linspace(0,Toptload,50);
            % %             sq.find('Keopsys MO').set(3.9); % do not touch the 3.9V of the Master oscillator or the conversion watts to V is ruined
            % %             sq.find('Keopsys FA').after(t,sq.minjerk(t,0,DipolePtoV('Keopsys',9))); %5
            % %             sq.find('Redpower TTL').set(1);
            % %             sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower',12))); %15
            %             sq.find('MOT bias').set(1);
            %             sq.find('MOT bias coil').after(t,sq.linramp(t,0,dipole_field));
            %             sq.delay(max(Toptload,Tmagload));
            t =  0:5e-3:Tmagload;
            dBLoad = 110; %80 %110 Ryan has that
            % sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',dBLoad/2),dBtoV('normal',dBLoad)));
            % sq.find('CD Fine/Fast').set(dBtoV('fine',10)); %Ryan Has 0
            sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',55),dBtoV('normal',dBLoad)));
            sq.find('CD Fine/Fast').set(dBtoV('fine',4));
            %    sq.find('CD Fine/Fast').set(6);
            %     sq.find('CD2').set(6); %6
            sq.delay(Tmagload); %Tmagload
            Toptload = 100e-3; %20e-3 %400e-3
            t = linspace(0,Toptload,50);
            sq.find('Keopsys FA').after(t,sq.minjerk(t,0,DipolePtoV('keopsys',8))); %8  %2.6 (when no PtoV used) %%%GOES TO 3.05 AT MOST % 7.5 W. 
            sq.find('Redpower TTL').set(1);
            sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower',12))); %14 % 15
            sq.find('MOT bias').set(1);
            sq.find('MOT bias coil').after(t,sq.linramp(t,0,dipole_field));
            sq.delay(max(Toptload,Tmagload));
%% 
% 
% 
% *Delay in the magtrap to see the life time (for testing purposes only!)*
% 
% Notes:
%% 
% * Using sq.delay(1) holds atoms in the mag trap for 1s
% * Using sq.delay(opt.params) holds atoms in the mag trap for a variable time. 
% This means going into the callback "Callback_MeasureTemperature_fancy.mlx" and 
% changing the variable to opt.params
            %             sq.delay(0); % delay for n second(s).
            % sq.delay(opt.params); %0 to as long as it can holds
%% 
% 
% 
% *Back to contents*
%% RF Evaporation
            if RF_Evaporation == 1
%% 
% 
                rf_start = 16; %12 Values in MHz % 12 prev % 16 assumed to be opt
                rf_end = 2; %1.6 if dipoles are correctly positioned
                rf_rate = 4; %MHz/s %4.55 % 2.5 prev
                Tevap = (rf_start - rf_end)/rf_rate;
                rf_ramp_type = 'lin';
                rf_exp_time_constant = 2;
                t = linspace(0,Tevap,100);
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
% *Test for Magnetic trap and optical dipole alignment (to be turned on only 
% in testing procedures- check within few ms tof).*
%                 sq.find('CD0 Fast').set(0);
%                 sq.find('CD Fine/Fast').set(0);
%                 sq.delay(20e-3 - opt.tof);
%% 
% 
% 
% 
% 
% *Turn off dipoles before ramping down coils (to be uncommented only for testing 
% procedure of the atoms in the magtrap)*
                % sq.find('Redpower CW').set(0);
                % sq.find('Redpower TTL').after(100e-6,0);
                % sq.find('Keopsys FA').set(0);
                % sq.find('Keopsys MO').after(100e-6,0);
%% 
% 
%% Load into Optical Dipoles - Ramp down coils
% 
                if Load_Opt_Evap == 1
%% 
% Note: 
%% 
% * For slow ramp down use 0.9 (s).
% * For fast ramp down use 0.3 (s).
                    Trampcoils = 0.4; %0.25 %0.4 opt
                    dB_weak = 0;
                    % rf_final = 1;
                    t = linspace(0,Trampcoils,51);
                    sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dBtoV('normal',dB_weak)));
                    sq.find('CD Fine/Fast').after(t,sq.linramp(t,sq.find('CD Fine/Fast').values(end),0));
%                     sq.find('CD Fine/Fast').set(dBtoV('fine',0));
                    sq.delay(Trampcoils);
%% 
% 
% 
% *Ask Ryan what this does!*
                    % % for bias coils
                    % sq.find('85 repump amp').set(0); % 3DMOT bias
                    % sq.find('MOT bias').set(0); %switching imaging coils on
                    % sq.find('MOT bias coil').set(0); %imaging coils
%                     to hold in dipole trap
%                     sq.delay(2); %1
%% 
% 
% 
% *Back to contents*
%% Optical Evaporation
                    if Optical_Evaporation == 1
%% 
% 
                        % Tevap = 2.5;   %3
                        Tevap = 3; % 3 opt
                        
%% 
% 
                        t = linspace(0,Tevap,250);
%                         final_dipole.RP = opt.redpower;
%                         final_dipole.FA = opt.keopsys;
% sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP),0.9)); %0.7
%                         sq.find('Keopsys FA').after(t,sq.expramp(t,sq.find('Keopsys FA').values(end),final_dipole.FA,0.8)); %0.7
                        final_dipole.RP = opt.redpower; % 1.00 opt
                        final_dipole.FA = opt.keopsys;  % 1.95 opt
                        sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP), 0.7)); %0.7 opt
%                         sq.find('Keopsys FA').after(t,sq.expramp(t,sq.find('Keopsys FA').values(end),final_dipole.FA,0.6)); %0.9 %0.7 *** THIS WAS ANABLED BEFORE, BUT IT DOES NOT REFERENCE THE KEOPSYS LASER
                        sq.find('Keopsys FA').after(t,sq.expramp(t,sq.find('Keopsys FA').values(end),DipolePtoV('keopsys',final_dipole.FA), 0.7)); %0.7 opt  
                        sq.find('variable wave plate').set(-2.42);
                        sq.delay(Tevap);
%% 
% 
%% Just hold in ODT:
% if(0)
%     sq.delay(opt.params);
% end 
%% 
% 
%% Adiabatically ramp up the laser powers:
if(1)
    %%% Define a time vector for optical evap: 
    TRamp = 70e-3; % 40ms
    t = linspace(0,TRamp,100);
    %%% Define the ODT params from 'opt' 
    final_dipole_ramp.RP = opt.param2;
    final_dipole_ramp.FA = opt.param1;
    
    %%% Set ODT params:
    sq.find('RedPower CW').after(t,sq.linramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole_ramp.RP))); 
    sq.find('Keopsys FA' ).after(t,sq.linramp(t,sq.find('Keopsys FA' ).values(end),DipolePtoV('keopsys', final_dipole_ramp.FA))); 
    
        %%% Delay for duration of ramping:
    sq.delay(TRamp);
    sq.delay(opt.params);
end 
%% 
% 
% 
% 
                    end % Optical evaporation
                end % Load in optical dipoles
            end % RF evaporation
        end %Load into the Mag trap
    end % PGC
end %end CMOT
%% 
% 
%% Drop those atoms
timeAtDrop = sq.time;
sq.find('Earth Bias 1').set(0);
sq.find('Earth Bias 2').set(0);
sq.find('Earth Bias 3').set(0);
sq.find('2D MOT Coils').set(0);
sq.find('3DMOT').set(0);
sq.find('87 repump').set(0);
sq.find('87 repump amp').set(0);
sq.find('CD0 Fast').set(dBtoV('normal',0)); %REMEMBER TO TURN BACK TO 2 WHEN USING RAMAN
sq.find('CD2').set(0);
sq.find('CD Fine/Fast').set(0);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
sq.find('Redpower CW').set(0);
sq.find('Redpower TTL').after(100e-6,0);
sq.find('Keopsys FA').set(0);
sq.find('Keopsys MO').after(100e-6,0);
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));
%% 
% *Back to contents*
%% Information Displayed
if  Display_info == 1
%% 
% 
    if CMOT == 1
        if PGC == 1
            if Load_Mag_Evap == 1
                if RF_Evaporation == 1
                    if Optical_Evaporation == 1
                        cprintf('Keywords','Optical evaporation Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
                        opt.stage = 'Optical evaporation';
                        baseOpt = evalin('base', 'opt');
                        baseOpt.stage = opt.stage;
                        assignin('base', 'opt', baseOpt);
                    else
                        cprintf('Keywords','RF evaporation Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
                        opt.stage = 'RF evaporation';
                        baseOpt = evalin('base', 'opt');
                        baseOpt.stage = opt.stage;
                        assignin('base', 'opt', baseOpt);
                    end
                else
                    cprintf('Keywords','Load into Mag Trap Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
                    opt.stage = 'Load into Mag trap';
                    baseOpt = evalin('base', 'opt');
                    baseOpt.stage = opt.stage;
                    assignin('base', 'opt', baseOpt);
                end
            else
                cprintf('Keywords','PGC Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
                opt.stage = 'PGC';
                baseOpt = evalin('base', 'opt');
                baseOpt.stage = opt.stage;
                assignin('base', 'opt', baseOpt);
            end
        else
            cprintf('Keywords','CMOT Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
            opt.stage = 'CMOT';
            baseOpt = evalin('base', 'opt');
            baseOpt.stage = opt.stage;
            assignin('base', 'opt', baseOpt);
        end
    else
        cprintf('Keywords','MOT Stage |Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
        opt.stage = 'MOT';
%         assignin('base','opt.stage',opt.stage);
    end
%% 
% *Back to contents*
%% Imaging Sequence
    sq.anchor(timeAtDrop);
    sq.camDelay = timeAtDrop - 2;
    makeImagingSequence(sq,'tof',opt.tof,'pulse time',40e-6,'repump delay',100e-6,...
        'repump time',200e-6,'cam time',5e-6,'cycle time',100e-3,...
        'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
        'fibre switch delay',1e-3,'imaging_field',imaging_field,'image type','horizontal'); %pulse time use to be 100e-6
    % makeImagingSequence_4_Images(sq,'tof',F2_tof,'tof2',F1_minus_F2_tof,'pulse time',100e-6,'repump delay',100e-6,...
    %     'repump time',200e-6,'cam time',50e-6,'cycle time',100e-3,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
    %     'fibre switch delay',1e-3,'imaging_field',imaging_field,'image type','horizontal');
    % turn off the dipoles
    sq.find('Redpower CW').set(0);
    sq.find('Redpower TTL').after(100e-6,0);
    sq.find('Keopsys FA').set(0);
    sq.find('Keopsys MO').after(100e-6,0);
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
%% 
% 
end %function