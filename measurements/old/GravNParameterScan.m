function GravNParameterScan(r)

if r.isInit()
    %Initialize run
    r.data.motamp = 2.75:0.25:3.5;
    r.data.motfreq = 4.75:0.25:5.5;
    
    r.c.setup('var',r.data.motamp,r.data.motfreq);
    r.numRuns = r.c.total;
    
    r.makerCallback = @makeSequence;
elseif r.isSet()
    %Increment counter
    if r.currentRun ~= 1
        r.c.increment;
    end
    
    %Build/upload/run sequence
    r.make(8.5,35e-3,1.8,r.data.motamp(r.c(1)),r.data.motfreq(r.c(2)));
    r.upload;
%     r.data.sq(r.currentRun,1) = r.sq.data;
    %Print information about current run
    fprintf(1,'Run %d/%d, MOT Amp: %.2f V, Freq: %.2f V\n',r.currentRun,r.numRuns,...
        r.data.motamp(r.c(1)),r.data.motfreq(r.c(2)));

elseif r.isAnalyze()
%     nn = r.currentRun;
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.5); %Wait for other image analysis program to finish with files
    
    %Analyze image data from last image
    c = Abs_Analysis('last');
    r.data.files{i1,i2} = {c.raw.files(1).name,c.raw.files(2).name};
    r.data.N(i1,i2) = c.N;
    r.data.Nth(i1,i2) = c.N.*(1-c.becFrac);
    r.data.Nbec(i1,i2) = c.N.*c.becFrac;
    r.data.F(i1,i2) = c.becFrac;
    r.data.xw(i1,i2) = c.gaussWidth(1);
    r.data.x0(i1,i2) = c.pos(1);
    r.data.yw(i1,i2) = c.gaussWidth(2);
    r.data.y0(i1,i2) = c.pos(2);
    r.data.OD(i1,i2) = c.fitdata.params.gaussAmp(1);
    r.data.Tx(i1,i2) = c.T(1);
    r.data.Ty(i1,i2) = c.T(2);
    
    figure(12);%clf;
    subplot(1,2,1);
    plot(r.data.motamp(1:i1),r.data.N(1:i1,i2),'o-');
    
    if r.c(1) == r.c.imax(1)
        subplot(1,2,2);
        cla;
        s = {};
        for nn = 1:i2
            plot(r.data.motamp,r.data.N(:,nn),'o-');
            s{nn} = sprintf('Freq %.2f V',r.data.motfreq(nn));
            hold on;
        end
        hold off;
        legend(s);
    end


end