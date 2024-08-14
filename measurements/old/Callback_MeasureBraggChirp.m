function Callback_MeasureBraggChirp(r)

if r.isInit()
    
%     r.data.chirp = const.randomize(-15e3:1e3:15e3);
    r.data.chirp = const.randomize(15e3:1e3:45e3);
    r.data.T = [5e-3,10e-3,15e-3];
    r.c.setup('var',r.data.chirp,r.data.T);
elseif r.isSet()
    
    r.make(r.devices.opt,'bragg',{'chirp',r.devices.opt.bragg.chirp + r.data.chirp(r.c(1)),'T',r.data.T(r.c(2))});
    r.upload;
    fprintf(1,'Run %d/%d, Chirp = %.3f MHz/s, T = %.3f ms\n',r.c.now,r.c.total,...
        r.data.chirp(r.c(1))/1e6,r.data.T(r.c(2))*1e3);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis('last');
    if ~img.raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif r.c.now > 1 && strcmpi(r.data.files{r.c.now - 1}.name,img.raw.files.name)
        r.c.decrement;
        return
    end
    %
    % Store raw data
    %
%     r.data.img{i1,1} = img;
    r.data.files{r.c.now} = img.raw.files;
    %
    % Get processed data
    %
    r.data.N(i1,i2,:) = img.get('N');
    r.data.Nsum(i1,i2,:) = img.get('Nsum');
    r.data.R(i1,i2) = r.data.N(i1,i2,1)./sum(r.data.N(i1,i2,:));
    r.data.Rsum(i1,i2) = r.data.Nsum(i1,i2,1)./sum(r.data.Nsum(i1,i2,:));

%     [~,N,dout] = FMI_Analysis;
%     r.data.N(i1,i2,:) = [N.N1,N.N2];
%     r.data.R(i1,i2) = N.R;
%     r.data.d{i1,i2} = dout;
%     
    figure(98);
    subplot(1,2,1);
    cla;
    h = plot(r.data.chirp(1:i1)/1e3,r.data.R(1:i1,i2),'o');
    for nn = 1:numel(h)
        set(h(nn),'MarkerFaceColor',h(nn).Color);
    end
    plot_format('Modified chirp [kHz/s]','Population','',12);
    grid on;
    ylim([0,1]);

    if r.c.done(1)
        subplot(1,2,2);
%         nlf = nonlinfit(r.data.chirp/1e6,r.data.R(:,i2));
%         nlf.setFitFunc(@(A,w,x0,x) A*exp(-(x - x0).^2/w^2));
%         nlf.bounds2('A',[0.4,0.7,0.5],'w',[0,1,0.05],'x0',[25.05,25.15,25.1]);
%         r.data.c{i2} = nlf.fit;
%         r.data.optimal_chirp(i2,1) = nlf.c(3,1);
        for nn = 1:i2
            h = plot(r.data.chirp/1e3,r.data.R(:,nn),'o');
            set(h,'MarkerFaceColor',h.Color);
            hold on
            grid on;
        end
        ylim([0,1]);
    end
    
end