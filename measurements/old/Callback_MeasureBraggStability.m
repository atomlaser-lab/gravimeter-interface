function Callback_MeasureBraggStability(r)

if r.isInit()
    
    r.c.setup(Inf);
elseif r.isSet()
    
%     r.make(25.5,730e-3,1.48,r.data.power(r.c(1)),0,130e-3);

%     r.make(r.devices.opt.set('power',r.data.power(r.c(1))));
    r.upload;
    fprintf(1,'Run %d/%d\n',r.c.now,r.c.total);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1);
%     img = Abs_Analysis('last');
%     if ~img(1).raw.status.ok()
%         %
%         % Checks for an error in loading the files (caused by a missed
%         % image) and reruns the last sequence
%         %
%         r.c.decrement;
%         return;
%     end
%     %
%     % Store raw data
%     %
% %     r.data.c{i1,1} = c;
%     r.data.files{i1,1} = img.raw.files;
%     %
%     % Get processed data
%     %
%     r.data.N(i1,:) = img.get('N');
%     r.data.Nsum(i1,:) = img.get('Nsum');
%     r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
%     r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));

    [~,N,dout] = FMI_Analysis;
    r.data.N(i1,:) = [N.N1,N.N2];
    r.data.R(i1,1) = N.R;
    r.data.Nsum(i1,:) = N.sum;
    r.data.Rsum(i1,1) = N.Rsum;
    r.data.dout(i1,1) = dout;
    
    r.devices.rp.getRAM;
    r.data.pulse(i1,1) = analyze_bragg_pulses(r.devices.rp.t,r.devices.rp.data(:,1),struct('width',24e-6));
    r.data.pulse(i1,2) = analyze_bragg_pulses(r.devices.rp.t,r.devices.rp.data(:,2),struct('width',24e-6));
    
    figure(98);clf;
    subplot(1,2,1);
    plot(1:i1,r.data.R(1:i1,:),'o-');
    hold on
    plot(1:i1,r.data.Rsum(1:i1,:),'sq-');
    plot_format('Run','Population','',12);
    grid on;
    ylim([0,1]);
    
    subplot(1,2,2);
    plot(1:i1,[r.data.pulse(:,1).pulse_area],'o-');
    hold on
    plot(1:i1,[r.data.pulse(:,2).pulse_area],'sq-');
    plot_format('Run','Pulse area','',12);
    grid on;
    
    
end