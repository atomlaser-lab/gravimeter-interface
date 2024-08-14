function GravMagLoad(r)

if r.isInit()
    %Initialize run
    r.data.param = const.randomize(1e-3:1e-3:10e-3);
    r.numRuns = numel(r.data.param);
    r.data.matlabfiles.callback = fileread('GravMagLoad.m');
    r.data.matlabfiles.init = fileread('gravimeter-interface/initSequence.m');
    r.data.matlabfiles.sequence = fileread('gravimeter-interface/makeSequenceFull.m');
    r.data.matlabfiles.analysis = fileread('Abs_Analysis.m');
elseif r.isSet()
    %Build/upload/run sequence
    r.make(r.data.param(r.currentRun),10e-3);
    r.upload;
    r.data.sq(r.currentRun,1) = r.sq.data;
    %Print information about current run
    fprintf(1,'Run %d/%d, Parameter: %.3f ms\n',r.currentRun,r.numRuns,...
        r.data.param(r.currentRun)*1e3);
    
elseif r.isAnalyze()
    nn = r.currentRun;
    pause(2.0); %Wait for other image analysis program to finish with files
    
    %Analyze image data from last image
    c = Abs_Analysis('last');
    r.data.files{nn,1} = c.raw.files(1).name;r.data.files{nn,2} = c.raw.files(2).name;
    r.data.N(nn,1) = c.N;
    r.data.xw(nn,1) = c.gaussWidth(1);
    r.data.x0(nn,1) = c.pos(1);
    r.data.yw(nn,1) = c.gaussWidth(2);
    r.data.y0(nn,1) = c.pos(2);
    r.data.OD(nn,1) = c.fitdata.params.gaussAmp(1);
    r.data.T(nn,:) = c.T;
    
    figure(10);clf;
    subplot(1,2,1);
    plot(r.data.param(1:nn),r.data.N/1e9,'o');
    plot_format('Bias voltage [V]','Number of atoms \times 10^9','',12);
    ylim([0,Inf]);
    grid on;
    subplot(1,2,2);
    plot(r.data.param(1:nn),r.data.OD,'o');
    plot_format('Bias voltage [V]','Peak OD','',12);
    ylim([0,Inf]);
    grid on;
end