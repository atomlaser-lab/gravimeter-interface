function Callback_MeasureRamanRamseyInterferometer(r)

p = 'time';

if r.isInit()
%     r.data.freq = -4.1 + 0.5 + const.randomize(500e-3*linspace(-1,1,21));
    r.data.(p) = 0.25:0.25:10;
    
    r.c.setup('var',r.data.(p));
elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.(p)(r.c(1)));
    r.upload;
    if strcmpi(p,'time')
        fprintf(1,'Run %d/%d, T = %.3f ms\n',r.c.now,r.c.total,r.data.(p)(r.c(1)));
    elseif strcmpi(p,'freq')
        fprintf(1,'Run %d/%d, F = %.3f kHz\n',r.c.now,r.c.total,r.data.(p)(r.c(1)));
    end
    
elseif r.isAnalyze()
    i1 = r.c(1);
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
    r.data.N(i1,:) = [img(1).get('N'),img(2).get('N')];
    r.data.Nsum(i1,:) = [img(1).get('Nsum'),img(2).get('Nsum')];
%     if numel(img) > 1
%         r.data.N(i1,2) = r.data.N(i1,2) - r.data.N(i1,1);
%         r.data.Nsum(i1,2) = r.data.Nsum(i1,2) - r.data.Nsum(i1,1);
%     end
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));
%     r.data.R2(i1,:) = r.data.N(i1,[3,7])./sum(r.data.N(i1,[3,7]));
    r.data.R2 = r.data.R;
    
    figure(97);clf;
    subplot(1,2,1)
    scatter(r.data.(p)(1:i1),r.data.N(1:i1,:),'filled');
    plot_format('Freq [kHz]','Number','',12);
    if strcmpi(p,'time')
        xlabel('Time [ms]');
    end
    title('Raman frequency using fit over OD')
    grid on
    hold on;
    ylim([0,Inf]);
    
    subplot(1,2,2)
    ax = scatter(r.data.(p)(1:i1),r.data.R2(1:i1,:),'filled');
    hold off;
    plot_format('Freq [kHz]','Population','',12);
    if strcmpi(p,'time')
        xlabel('Time [ms]');
    end

    grid on;
%     if r.c.done
%         tNow = datestr(now);
%         caption = sprintf('Determination of Raman frequency %s', tNow);
%         sgtitle(caption)
%     end
    if r.c(1) > 5
        nlf = nonlinfit(r.data.(p)(1:r.c(1)),r.data.R2(:,2),1e-2);
        if strcmpi(p,'freq')
            nlf.setFitFunc(@(y0,C,f,x) y0 + C/2*cos(2*pi*(x - f)*1));
        elseif strcmpi(p,'time')
            nlf.setFitFunc(@(y0,C,f,ph,x) y0 + C/2*cos(2*pi*f*x + ph));
        end
        [~,idx] = min(nlf.y);
        nlf.bounds2('y0',[0,1,0.5],'C',[0,1,0.8],'f',[0,1,.1],'ph',[-2*pi,2*pi,0]);
        nlf.fit;
        subplot(1,2,2);
        hold on
        xplot = linspace(min(nlf.x),max(nlf.x),1e2);
        plot(xplot,nlf.f(xplot),'-');
        hold off;
        if r.c.done
            nlf.montecarlo(100,1);
            fprintf(1,'Detuning: %.3f +/- %.3f kHz\n',nlf.c(3,:));
        else
            fprintf(1,'Detuning: %.3f kHz\n',nlf.c(3,1));
        end
    end
end