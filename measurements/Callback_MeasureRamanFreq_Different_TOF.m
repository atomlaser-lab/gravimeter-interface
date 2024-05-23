function Callback_MeasureRamanFreq_Different_TOF(r)

if r.isInit()
    r.data.pulse_duration = 5;
    r.data.freq = 0e-3 + const.randomize(250*linspace(-1,1,15));
    r.data.tof = [0.1e-3,0.5e-3,1e-3,2e-3,3e-3,5e-3,10e-3];
    
    r.c.setup('var',r.data.freq,r.data.tof);
elseif r.isSet()
    r.make(r.devices.opt,'params',[r.data.pulse_duration,r.data.freq(r.c(1)),r.data.tof(r.c(2))]);
    r.upload;
    fprintf(1,'Run %d/%d, F = %.3f kHz, TOF = %.1f ms\n',r.c.now,r.c.total,...
        r.data.freq(r.c(1)),r.data.tof(r.c(2))*1e3);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.5 + 0.5*rand);
    img = Abs_Analysis_DualState_RT('last');
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    end
    %
    % Store raw data
    %
    r.data.files{i1,1} = img(1).raw.files;
    %
    % Get processed data
    %
    r.data.N{i2}(i1,:) = [img(1).get('N'),img(2).get('N')];
    r.data.Nsum{i2}(i1,:) = [img(1).get('Nsum'),img(2).get('Nsum')];

    r.data.R{i2}(i1,:) = r.data.N{i2}(i1,:)./sum(r.data.N{i2}(i1,:));
    r.data.Rsum{i2}(i1,:) = r.data.Nsum{i2}(i1,:)./sum(r.data.Nsum{i2}(i1,:));
    r.data.R2{i2}(i1,:) = r.data.N{i2}(i1,[3,7])./sum(r.data.N{i2}(i1,[3,7]));
    
    figure(97);
    subplot(1,2,1);cla;
    scatter(r.data.freq(1:i1),r.data.R2{i2}(1:i1,:),'filled');
    plot_format('Freq [kHz]','Number','',12);
    grid on
    hold on;
    ylim([0,Inf]);

    if r.c.done(1)
        subplot(1,2,2);
        scatter(r.data.freq(1:i1),r.data.R2{i2}(1:i1,:),'filled');
        hold on;
        plot_format('Freq [kHz]','Population','',12);
        grid on;
  
        nlf = nonlinfit(r.data.freq(1:r.c(1)),r.data.R2{i2}(:,2),1e-2);
        nlf.setFitFunc(@(A,R,x0,x) A*(1 - 4*R.^2./(4*R^2+(x-x0).^2).*sin(2*pi*sqrt(4*R^2+(x-x0).^2).*5e-3/2).^2));
        [~,idx] = min(nlf.y);
        nlf.bounds2('A',[0.5,2,0.95],'R',[0,30,0.3],'x0',[min(nlf.x),max(nlf.x),nlf.x(idx)]);
        nlf.ex = sum(r.data.N{i2},2) < 0.5*max(sum(r.data.N{i2},2));
        nlf.fit;
        fprintf(1,'Rabi frequency: %.3f kHz, Center = %.3f kHz\n',nlf.c(2,1),nlf.c(3,1));
        r.data.rabi(i2,1) = nlf.c(2,1);
        r.data.center_freq(i2,1) = nlf.c(3,1);

        xplot = linspace(min(nlf.x),max(nlf.x),1e2);
        plot(xplot,nlf.f(xplot),'-','HandleVisibility','off');

        str = strsplit(strtrim(sprintf('%.1f ms\n',1e3*r.data.tof(1:i2))),'\n');
        legend(str);
    end
end