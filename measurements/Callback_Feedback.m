function Callback_Feedback(r)

if r.isInit()
%     r.data.param = logical(mod(1:1000,2));
    r.data.param = const.randomize(20:50);
    r.c.setup('var',r.data.param);
elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.param(r.c(1))).upload;
    fprintf(1,'Run %d/%d, Param = %.3f\n',r.c.now,r.c.total,r.data.param(r.c(1)));
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.25*rand);
    [img,nd,fb] = Abs_Analysis('last',1);
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end
    
    r.data.files{i1,1} = img(1).raw.files;
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');
    r.data.becFrac(i1,:) = img.get('becFrac');
    r.data.OD(i1,:) = img.get('peakOD');
%     r.data.T(i1,:) = img.get('T');
    Ttmp = img.get('T');
    r.data.T1(i1,:) = Ttmp(:,:,1);
    r.data.T2(i1,:) = Ttmp(:,:,2);
%     r.data.pos(i1,:) = img.get('pos');
    
    r.data.xwidth(i1,1) = var(fb.xwidth(1:75));
    r.data.xwidth(i1,2) = var(fb.xwidth(76:150));
    r.data.fb(i1) = fb;


%     figure(98);clf;
%     subplot(1,2,1);
%     plot(r.data.param(1:i1),r.data.xwidth,'o');
%     grid on
%     plot_format('Driving Frequency [Hz]','X variance [pixels^2]','',10);
%     subplot(1,2,2);
%     plot(r.data.param(1:i1),r.data.T*1e9,'o');
%     grid on;
%     plot_format('Driving Frequency [Hz]','Temperature [nK]','',10);
%     legend('T_x','T_y');

    figure(98);clf;
    subplot(1,2,1);
    plot(r.data.param(1:i1),r.data.xwidth,'o');
    grid on
    plot_format('Run','X variance [pixels^2]','',10);
    subplot(1,2,2);
    plot(r.data.param(1:i1),r.data.T2*1e9,'o');
    grid on;
    plot_format('Run','Temperature [nK]','',10);
    legend('T_x','T_y');

%     subplot(1,2,1);
%     plot(1:i1,r.data.T1*1e9,'o');
%     grid on
%     plot_format('Run','Temperature [nK]','|1,-1> temperature',10);
%     ylim([0,200]);
%     subplot(1,2,2);
%     plot(1:i1,r.data.T2*1e9,'o');
%     grid on
%     plot_format('Run','Temperature [nK]','|1,0> temperature',10);
%     ylim([0,200]);

%     figure(99);clf;
%     idx = r.data.param(1:i1);
%     subplot(2,2,1);
%     histogram(r.data.T1(idx,1)*1e9,60:200)
%     hold on
%     histogram(r.data.T1(~idx,1)*1e9,60:200)
%     grid on
%     plot_format('T_x [nK]','Counts','|1,-1> T_x',10);
%     legend('Mod on','Mod off')
%     legend({'Mod on','Mod off'},'location','northwest')
%     subplot(2,2,2);
%     histogram(r.data.T1(idx,2)*1e9,60:200)
%     hold on
%     histogram(r.data.T1(~idx,2)*1e9,60:200)
%     grid on
%     plot_format('T_y [nK]','Counts','|1,-1> T_y',10);
%     subplot(2,2,3);
%     histogram(r.data.T2(idx,1)*1e9,60:200)
%     hold on
%     histogram(r.data.T2(~idx,1)*1e9,60:200)
%     grid on
%     plot_format('T_x [nK]','Counts','|1,0> T_x',10);
%     subplot(2,2,4);
%     histogram(r.data.T2(idx,2)*1e9,60:200)
%     hold on
%     histogram(r.data.T2(~idx,2)*1e9,60:200)
%     grid on
%     plot_format('T_y [nK]','Counts','|1,0> T_y',10);

end


end