% % % Inputs

fid = fopen('C:\Users\Apollo\Desktop\Temp\PhaseMeasurement_04072024\RyanSystem1.BIN','rb'); 

Title = 'Phase Noise With Power Lock';

dt = 1/1.09e6; %1/sampleRate of Scope 
freq = 2*pi*((112.50000047-110.000000102)*2)*1e6; %2*pi*Expected beat frequency
FigNum = 5;

%% Load

fseek(fid,0,'eof');
fsize = ftell(fid);
frewind(fid);
raw = fread(fid,fsize,'uint8');
fclose(fid);

%% Process
data = zeros(1e6,1,'int16');
mm = 1;
for nn = 1:2:(size(data,1)*2)
%     rr = typecast(uint8(raw(nn + [1,0])),'int8');
%     rr = rr/2;
%     rr = typecast(rr,'uint8');
%     data(mm) = typecast(rr,'int16');
    data(mm) = typecast(uint8(raw(nn + [1,0])),'int16');
    
    mm = mm + 1;
end
data = double(data);
data(data < 0) = data(data < 0) - 2*double(intmin('int16'));
Nsamples = numel(data);

% dt = 1/31.25e6;
t = dt*(0:1:(Nsamples - 1))';

%% Get phase


Iraw = data.*sin(freq*t);
Qraw = data.*cos(freq*t);

R = 2^6;
[I,tmeas] = cicfilter(t,Iraw,R,3);
[Q,~] = cicfilter(t,Qraw,R,3);

ph = unwrap(atan2(Q,I));

%% Plot
figure(FigNum);clf;
subplot(2,1,1);hold on
plot(tmeas,ph,'.-')
plot_format('Time [s]','Phase [rad]','',10);
grid on;
sgtitle(Title)

subplot(2,1,2);hold on
const.plotfft(detrend(ph,1),tmeas,[],'psd');

function [x,t] = cicfilter(t,x,R,N)

%% Integrator
xx = x;
for nn = 1:N
    xx = cumsum(xx,1);
end
%% Rate reduction
tr = t(R:R:end);
xr = xx(R:R:end,:);
%% Comb
xx = xr;
for nn = 1:N
    xx = diff(xx,1,1);
end
t = tr(1:size(xx,1));
x = R.^-N.*xx;

end


