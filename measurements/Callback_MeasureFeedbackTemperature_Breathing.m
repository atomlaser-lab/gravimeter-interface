function Callback_MeasureFeedbackTemperature_Breathing(r)

if r.isInit()
    r.data.runs = 1 : 50;
    r.data.enable = logical(mod(r.data.runs,2));
%     r.data.delay = 30e-3*rand(numel(r.data.runs),1);

    r.c.setup('var',r.data.runs);

    r.devices.fb = FeedbackControl;
    r.devices.fb.gains = [0,0,0,0;0,-0.5,0,0]; % WAS [-0.5,0,0,0;0,-0.4,0,0];
    r.devices.fb.start = 20;
    r.devices.fb.stop = 100;
    r.devices.fb.jump_enable = 0;
    r.devices.fb.jumps = [0;0;0;0];
    r.devices.fb.jump_index = 10;
    r.devices.fb.biases = [6.3;0];
elseif r.isSet()
    r.devices.fb.enable = r.data.enable(r.c(1));
    % r.devices.fb.enable = false;
    r.devices.fb.upload;
    r.devices.opt.detuning = 6;
    r.devices.opt.load_time = 15;
    r.devices.opt.raycus = 1.0;
    r.devices.opt.redpower = 1.1;
    r.make(r.devices.opt,'params',r.data.enable(r.c(1))).upload;
    % fprintf(1,'Run %d/%d\n',r.c.now,r.c.total);
    % fprintf('Enable = %d\n',r.data.enable(r.c(1)));
    if r.data.enable(r.c(1)) == 1
      fprintf(1,'Run %d/%d FB enabled',r.c.now,r.c.total);
    elseif r.data.enable(r.c(1)) == 0
      fprintf(1,'Run %d/%d FB disabled',r.c.now,r.c.total);  
    end 
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.25*rand);
    [img, nd, fbc] = Abs_Analysis_FB('last',1);
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(img.raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        pause(15);
        return;
    end
    
    r.data.files{i1,1} = img.raw.files;
    r.data.N(i1,1) = img.get('N');
    r.data.becFrac(i1,1) = img.get('becFrac');
    r.data.OD(i1,1) = img.get('peakOD');
    r.data.T(i1,1) = prod(squeeze(img.get('T')))^0.5;
    r.data.x(i1,1) = img.clouds.pos(1);
    r.data.y(i1,1) = img.clouds.pos(2);

    %%% In trap
    r.data.x_inTrap(i1,:)  = fbc.xpos(:);
    r.data.z_inTrap(i1,:)  = fbc.ypos(:);
    r.data.wx_inTrap(i1,:) = fbc.xwidth(:);
    r.data.wz_inTrap(i1,:) = fbc.ywidth(:);

    %%% Store feedback structure:
    r.data.fbc(i1) = fbc;


%     data = r.data;save('D:\data\22-11-19\collected-data-19-11-2022','data');

    if i1 > 5 && all(r.data.N((i1-5):i1) < 0.8e5)
        r.c.final = i1;
    end

    xx = r.data.runs(1:i1);
    x = xx(r.data.enable(1:i1));
    y = xx(~r.data.enable(1:i1));
    figure(98);clf;
    subplot(1,3,1);
    h = plot(x,r.data.N(r.data.enable(1:i1)),'o');
    set(h,'MarkerFaceColor',h.Color);
    hold on
    h = plot(y,r.data.N(~r.data.enable(1:i1)),'sq');
    set(h,'MarkerFaceColor',h.Color);
    ylim([0,Inf]);
    grid on
    plot_format('Run','Number','',10);
    legend('FB On','FB Off')

    subplot(1,3,2);
    h = plot(x,r.data.T(r.data.enable(1:i1))*1e9,'o');
    set(h,'MarkerFaceColor',h.Color);
    hold on
    h = plot(y,r.data.T(~r.data.enable(1:i1))*1e9,'sq');
    set(h,'MarkerFaceColor',h.Color);
    % ylim([80,200]);
    grid on
    plot_format('Run','Temperature [nK]','',10);
    legend(sprintf('FB On %.0f +/- %.0f nK', 1e9 * mean(r.data.T(r.data.enable(1:i1))), 1e9 * std(r.data.T(r.data.enable(1:i1))) / sqrt(length(1:i1)/2)), sprintf('FB Off %.0f +/- %.0f nK', 1e9 * mean(r.data.T(~r.data.enable(1:i1))), 1e9 * std(r.data.T(~r.data.enable(1:i1))) / sqrt(length(1:i1)/2)));

    subplot(1,3,3);
    h = plot(x,r.data.becFrac(r.data.enable(1:i1))*1e2,'o');
    set(h,'MarkerFaceColor',h.Color);
    hold on
    h = plot(y,r.data.becFrac(~r.data.enable(1:i1))*1e2,'sq');
    set(h,'MarkerFaceColor',h.Color);
    % ylim([80,200]);
    grid on
    plot_format('Run','BEC fraction [percent]','',10);
    legend(sprintf('FB On %.0f +/- %.0f per', 1e2 * mean(r.data.becFrac(r.data.enable(1:i1))), 1e2 * std(r.data.becFrac(r.data.enable(1:i1))) / sqrt(length(1:i1)/2)), sprintf('FB Off %.0f +/- %.0f per', 1e2 * mean(r.data.becFrac(~r.data.enable(1:i1))), 1e2 * std(r.data.becFrac(~r.data.enable(1:i1))) / sqrt(length(1:i1)/2)));
end


end