function GravInterferometerTest(r)

if r.isInit()
%     r.data.chirp0 = 2*r.devices.p.k*r.devices.p.g/(2*pi);
%     r.data.param = r.data.chirp0 + (-100e3:10e3:100e3);
    r.data.param = 0:0.005:0.15;

    r.numRuns = numel(r.data.param);
elseif r.isSet()
%     r.devices.p.makeSinglePulse(r.data.param(r.currentRun));

%     r.devices.p.makePulses('width',r.data.param(r.currentRun),'braggpower',[r.data.parame(r.currentRun),0,0]);
%     r.devices.p.makePulses('chirp',r.data.param(r.currentRun));
    r.devices.p.makePulses('power',[r.data.param(r.currentRun),0,0]);

%     tic;
%     r.devices.p.upload;
%     toc;
    r.make(8.25,35e-3,1.1);
    r.upload;
%     fprintf(1,'Run %d/%d, Chirp: %.2f\n',r.currentRun,r.numRuns,r.data.param(r.currentRun) - r.data.chirp0);
    fprintf(1,'Run %d/%d, Power: %.2f\n',r.currentRun,r.numRuns,r.data.param(r.currentRun));
    
elseif r.isAnalyze()
    nn = r.currentRun;
    pause(1);
    c = Abs_Analysis_NClouds('last');
%     v = Abs_Analysis('last');
%     r.data.files{nn,1} = c.raw.files(1).name;r.data.files{nn,2} = c.raw.files(2).name;
    r.data.c(nn,:) = c;
    r.data.N(nn,:) = [c.N];
    r.data.R(nn,1) = r.data.N(nn,1)./sum(r.data.N(nn,:));
%     r.data.v(nn,:) = v;
    
    figure(10);clf;
%     chirp = r.data.param(1:nn)-r.data.chirp0;
    subplot(1,2,1);
    for mm = 1:size(r.data.N,2)
        errorbar(r.data.param(1:nn),r.data.N(1:nn,mm)/1e6,0.05*r.data.N(1:nn,mm)/1e6+0.05,'o');
        hold on;
    end
    subplot(1,2,2);
    errorbar(r.data.param(1:nn),r.data.R,0.02*ones(size(r.data.R)),'o');
    pause(0.01);
%     figure(11);clf;
%     for j=1:nn
%         subplot(8,10,j)
%         r.data.v(j).plotAbsData([0,.1],1);
% %         title(sprintf('chirp: %.3f',r.data.param(j)*1e6));
%     end
%     pause(0.01);
    
end