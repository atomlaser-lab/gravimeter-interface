clear
clc

%NOTE: make sure that the "last-image.txt" file has the corresponding last
%image number

% addpath("cloud images\");
% addpath("imaging-analysis\");

number_of_images = 1950; %141

% last_image = load('last-image.txt');
counter = 0;

for jj = 1:number_of_images

    %     FileInfo = dir(sprintf('C:\\Users\\QSensors\\OneDrive - Australian National University\\Documents\\GingerData\\cloud images\\bec%d.bin',last_image - jj + 1));
    %     data.time(jj,:) = FileInfo.date;
    % [data.year(jj,:),data.month(jj,:),data.day(jj,:),data.hour(jj,:),data.minute(jj,:),data.second(jj,:)] = datevec(FileInfo.datenum);
    % img = Abs_Analysis_DualState("last",jj);
    img = Abs_Analysis_DualState("last",174046 - 170271 - jj);


    data.files{jj,1} = img(1).raw.files;
    data.N(jj,:) = img.get('N');
    data.Nsum(jj,:) = img.get('Nsum');
    data.peakOD(jj,:) = img.get('peakOD');
    data.R(jj,:) = data.N(jj,:)./sum(data.N(jj,:));
    data.Rsum(jj,:) = data.Nsum(jj,:)./sum(data.Nsum(jj,:));

    if (~img(1).raw.status.ok() || data.N(jj,1) == 0 || data.N(jj,1) > 1e7 || data.Rsum(jj,1) - data.Rsum(jj,2) > 0.9 || data.Rsum(jj,1) - data.Rsum(jj,2) < 0.1 || isnan(data.Rsum(jj,1)))
    a = 1;
    else
        counter = counter + 1;
        Rsum(counter,:) = data.Rsum(jj,:);

        % plotting interferometer data
        figure(3)
        clf
        scatter([1:counter],Rsum(:,1) - Rsum(:,2),'filled')
        axis square
        xlabel('Number of Runs')
        ylabel('N_1-N_2')
        title('Allan Deviation Data')

        % Allan deviation

        raw_pop_data = Rsum; % insert your data here
        N1 = raw_pop_data(:,1);
        N2 = raw_pop_data(:,2);

        contrast =  0.1207*3; %0.07063
        DC_offset = 0.4122; %0.7318
        phase_offset = 9.51;
        l = 10;
        T = 1e-3;

        duty_cycle = 16; %this is the duty cycle of your machine
        Fs = 1/duty_cycle;
%         cosphi = asin((N1-N2-DC_offset)/(2*contrast)) + phase_offset;
cosphi = asin((1+N1-N2-2*DC_offset)/contrast) + phase_offset;
        % cosphi = asin((N1-N2-offset)/contrast)/(2*1e-3*10);
        %Convert the population to a phase
        phase = cosphi;
        phase = phase(~isnan(phase));

        if length(phase) > 2 && imag(cosphi(counter)) == 0

            %Averaging factod
            % maxFactor = floor(length(phase) / 2); % Maximum averaging factor
            % m = unique(ceil(logspace(0, log10(maxFactor), 20))); %averaging factor to work with
            maxFactor = floor(length(phase) / 2); % Maximum averaging factor
            % m = unique(ceil(logspace(0, log 10(maxFactor), 20)));
            m = 'octave';

            % Calculate Allan variance for each data set
            [av_data, tau_data] = allanvar(phase,m,Fs);

            % Introduce the averaging factor
            m_estimated = floor(length(phase) ./ (Fs * tau_data)); % u might want to adjust to get vector column for it to work

            % Calculate Allan deviation (square root of allan variance)
            ad_data= sqrt(av_data);

            % Calculate the proportional 1/√N trend
            tau_trend = tau_data; % Use the tau values from the Adev calculation
            trend_line = ad_data(1)*sqrt(duty_cycle)./ sqrt(tau_trend);

            % Calculate the error bars for the Adev
            error_data = ad_data ./ sqrt(m_estimated); %maybe should be

            % Plot Allan deviation
            figure(991)
            clf;
            %This returns the data points of the allan dev wiht errorbars an in runs (not in time -remove /duty cylce for itime)

            errorbar(tau_data/duty_cycle, ad_data*1000, error_data, '*-');
            hold on;
            loglog(tau_data/duty_cycle, trend_line*1000, 'k--');
            xlabel('Runs')
            ylabel('Allan Deviation ( mrad)')

            % Set the axes to logarithmic scale
            ax = gca;
            set(ax, 'XScale', 'log', 'YScale', 'log')
            grid off;
            title('Allan Deviation Phase')
        end
    end
end

% figure(992)
% clf;
% %This returns the data points of the allan dev wiht errorbars an in runs (not in time -remove /duty cylce for itime)
%
% errorbar(tau_data/duty_cycle, ad_data*(1/(2*l*T)), error_data, '*-');
% hold on;
% loglog(tau_data/duty_cycle, trend_line*(1/(2*l*T)), 'k--');
% xlabel('Runs')
% ylabel('Allan Deviation ( rad)')
% title('Allan Deviation Rotation Rate')
%
% % Set the axes to logarithmic scale
% ax = gca;
% set(ax, 'XScale', 'log', 'YScale', 'log')
% grid off;
