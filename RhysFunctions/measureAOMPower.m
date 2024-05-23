function data = measureAOMPower(ch)
%
% Short script for calibrating the Bragg AOMS.  Takes a mogchannel object
% and writes different power values to that.  The user enters the measured
% powers on the command prompt.
%
% Use: 
% 1) create mogdevice object: mog = mogdevice;mog.connect('192.168.1.4');
% 2) create 2 mogchannel objects ch1 = mogchannel(mog,1);ch2 = mogchannel(mog,2);
% 3) set to standard mode using: ch1.read;ch2.read;
% 4) Run for each channel, e.g.: data1 = measureAOMPower(ch1);
%
if strcmpi(ch.powunits,'dbm')
    P_applied = [-50,0,2:2:16]';
    Prf = 1e-3*10.^(P_applied/10);
    Popt = zeros(numel(Prf),1);
elseif strcmpi(ch.powunits,'hex')
    P_applied = [0:500:10000,11000:1000:16000]';
    Popt = zeros(numel(P_applied),1);
else
    error('Power units for the channel can only be ''dBm'' or ''hex''');
end

for nn = 1:numel(Popt)
    try
        ch.write('signal',1,'amplifier',1,'power',P_applied(nn));
    catch err
        return;
    end
    if strcmpi(ch.powunits,'dbm')
        fprintf(1,'RF Power = %.3f dBm\n',P_applied(nn));
    else
        fprintf(1,'RF amplitude = %d\n',round(P_applied(nn)));
    end
    Popt(nn,1) = input('Measured Power: ');
end

ch.write('amplifier',0,'power',0);

if strcmpi(ch.powunits,'dbm')
    data.dB = P_applied;
    data.Prf = Prf;
    data.Popt = Popt;
    x = Prf;
    nlf = nonlinfit(x,data.Popt - data.Popt(1),20e-6+1e-2*data.Popt);
    nlf.setFitFunc(@(y0,A,x0,x)y0 + A*sin((x./x0).^0.5.*pi/2).^4);
    nlf.bounds([-50e-6,0,0],[100e-6,1,4],[0,0.4,2]);
    figure(5);clf;
    nlf.fit,nlf.plot;
    data.nlf = nlf;
    data.rfscale = nlf.c(3,1);
else
    data.amp = P_applied;
    data.optical_power = Popt;
    x = P_applied;
end
