function GravTempMeas(r)

if r.isInit()
    %Initialize run
    r.data.tof = 25e-3:5e-3:35e-3; 
    r.data.param= 1:50; 
    r.c.setup('var',r.data.tof,r.data.param);
    
elseif r.isSet()
    r.make('detuning',0,'tof',r.data.tof(r.c(1)),'dipole',2.25,'power',0.,'T',10e-3,'camera','drop 1');

    r.upload;

    %Print information about current run
    fprintf(1,'Run %d/%d, TOF: %.1f ms, Param: %.3f\n',r.c.now,r.c.total,...
        r.data.tof(r.c(1))*1e3,r.data.param(r.c(2)));
    
elseif r.isAnalyze()
  

    i1=r.c(1);
    i2=r.c(2);
    
    pause(0.5);
    
    c = Abs_Analysis('last',1,'tof',r.data.tof(r.c(1)));
%     r.data.N(i1,i2,:) = reshape(c.get('N'),[1,1,2]);
%     r.data.Nsum(i1,i2,:) = reshape(c.get('Nsum'),[1,1,2]);
    r.data.N(i1,i2) = c.get('N');
    r.data.T(i1,i2,:) = reshape(c.clouds.T,[1,1,2]);
    r.data.w(i1,i2,:) = reshape(c.clouds.gaussWidth,[1,1,2]);
    r.data.OD(i1,i2) = c.clouds.peakOD;

    
        if ~c(1).raw.status.ok() 
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    end
    figure(10);
    subplot(1,3,1);
    w = squeeze(r.data.w(1:i1,i2,:));
    plot(r.data.tof(1:i1),w,'o');
    plot_format('Param [V]','Widths [m]','',12);
    grid on;
    if r.c.done(1)
        subplot(1,3,2);
        cla;
        errorbar(r.data.param(1:i2),r.data.N(i1,1:i2),0.025*r.data.N(i1,1:i2),'o');
        
        figure(11);clf;
        for nn = 1:2
            w = squeeze(r.data.w(:,i2,nn));
            lf = linfit(r.data.tof,w.^2,2*w.*20e-6);
            lf.setFitFunc('poly',[0,2]);
            lf.fit;
            lf.plot;
            r.data.Tfit(i2,nn) = lf.c(2,1)*const.mRb/const.kb*1e6;
            r.data.Terr(i2,nn) = lf.c(2,2)*const.mRb/const.kb*1e6;
        end
        
        figure(10);
        subplot(1,3,3);
        cla;
        for nn = 1:2
            errorbar(r.data.param(1:i2)',r.data.Tfit(1:i2,nn),r.data.Terr(1:i2,nn),'o');
            hold on;
        end
        ylim([0,Inf]);
        fprintf(1,'Tx = %.3f, Ty = %.3f\n',r.data.Tfit(i2,:));
    end
end


end

