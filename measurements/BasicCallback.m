function BasicCallback(r)

if r.isInit()
    r.data.param = const.randomize(5:15);
    r.c.setup('var',r.data.param);
elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.param(r.c(1))).upload;
    fprintf(1,'Run %d/%d, Param = %.2f\n',r.c.now,r.c.total,r.data.param(r.c(1)));
elseif r.isAnalyze()
    i1 = r.c(1);
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

    r.data.files{i1} = img.raw.files;
    r.data.N(i1,1) = squeeze(img.get('N'));
    r.data.F(i1,1) = squeeze(img.get('becFrac'));
    r.data.T(i1,1) = sqrt(prod(img.get('T')));
    r.data.pos(i1,:) = img.clouds(1).pos;
    r.data.width(i1,:) = img.clouds(1).gaussWidth;
    figure(10);clf;
    plot(r.data.param(1:i1),r.data.N(1:i1,:),'o');
    ylim([0,Inf]);
    grid on;
end


end