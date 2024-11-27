function Callback_LoadRate(r)

if r.isInit()
    r.data.load_time = [0.25,0.5,1:10,20:20:60]; % Load time
    r.c.setup('var',r.data.load_time);
elseif r.isSet()
    r.make(r.devices.opt,'load_time',r.data.load_time(r.c(1))).upload;
%     r.make(r.devices.opt).upload;
    fprintf(1,'Run %d/%d, Load time = %.2f s\n',r.c.now,r.c.total,r.data.load_time(r.c(1)));
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.25*rand);
    img = Abs_Analysis('last',1);
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
%     r.data.T(i1,:) = prod(squeeze(img.get('T')))^0.5;
%     r.data.pos(i1,:) = img.get('pos');

    figure(98);clf;
    plot(r.data.load_time(1:i1),r.data.N,'o');
    ylim([0,Inf]);
    grid on
    plot_format('Loading time [s]','Atom number','',10);
    if i1 > 3
        nlf = nonlinfit(r.data.load_time(1:numel(r.data.N)),r.data.N);
        nlf.setFitFunc(@(A,D,x) A*(1 - exp(-D*x)));
        nlf.bounds2('A',[0,5*max(nlf.y),max(nlf.y)],'D',[0,5,0.01]);
%         nlf.ex = (nlf.y == 0) | (nlf.x > 8);
        nlf.fit;
        hold on
        plot(nlf.xplot,nlf.f(nlf.xplot),'--');
        text(0.5,0.5,sprintf('Load rate = %.1e atoms/s\nDecay time = %.2f s\nMax number = %.1e\n',nlf.c(1,1)*nlf.c(2,1),1/nlf.c(2,1),nlf.c(1,1)),'units','normalized');
        fprintf('Load rate = %.1e, Decay time = %.2f, Max number = %.1e\n',nlf.c(1,1)*nlf.c(2,1),1/nlf.c(2,1),nlf.c(1,1));
    end

end


end