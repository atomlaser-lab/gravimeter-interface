function Callback_MeasureG(r)

if r.isInit()
    r.data.phase = 0:20:180;
    r.data.T = [1,2,3,4,5]*1e-3;
    r.c.setup('var',r.data.phase,r.data.T);
elseif r.isSet()
    
    r.make(r.devices.opt,'bragg',{'phase',r.data.phase(r.c(1)),'T',r.data.T(r.c(2))});
    r.upload;
    fprintf(1,'Run %d/%d, T = %.2f ms, Phase: %.2f\n',r.c.now,r.c.total,...
        r.data.T(r.c(2))*1e3,r.data.phase(r.c(1)));
    
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
    
    r.data.files{r.c.now} = img.raw.files;
    %
    % Get processed data
    %
    r.data.N(i1,i2,:) = img.get('N');
    r.data.Nsum(i1,i2,:) = img.get('Nsum');
    r.data.R(i1,i2) = r.data.N(i1,i2,1)./sum(r.data.N(i1,i2,:));
    r.data.Rsum(i1,i2) = r.data.Nsum(i1,i2,1)./sum(r.data.Nsum(i1,i2,:));
    
    nlf = nonlinfit;
    nlf.setFitFunc(@(y0,A,phi,x) y0+A/2*cos(2*x+phi));
    nlf.bounds2('y0',[0,1,0.5],'A',[0,1,0.7],'phi',[-pi,pi,0]);
    
    figure(98);
    subplot(1,2,1);
    errorbar(r.data.phase(1:i1),r.data.R(1:i1,i2),0.02*ones(size(r.data.R(1:i1,i2))),'o');
    hold on;
    errorbar(r.data.phase(1:i1),r.data.Rsum(1:i1,i2),0.02*ones(size(r.data.Rsum(1:i1,i2))),'sq');
    hold off;
    plot_format('Phase [deg]','N_{rel}','',12);
    subplot(1,2,2);
    if r.c.done(1)
        cla;
        for mm = 1:i2
            errorbar(r.data.phase,r.data.Rsum(:,mm),0.02*ones(size(r.data.Rsum(:,mm))),'o');
            hold on;
            nlf.set(r.data.phase*pi/180,r.data.Rsum(:,mm),0.02*ones(size(r.data.Rsum(:,mm))));
            r.data.coeffs{mm} = nlf.fit;
            plot(r.data.phase,nlf.f(r.data.phase*pi/180),'-');
            str{mm} = sprintf('Time = %d ms',round(1e3*r.data.T(mm)));
        end
        plot_format('Phase [deg]','N_{rel}','',12);
        legend(str);
    end
    pause(0.01);
    
end