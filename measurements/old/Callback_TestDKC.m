function Callback_TestDKC(r)

if r.isInit()
    r.data.dt = const.randomize(1e-3:1e-3:10e-3);
    
    r.c.setup('var',r.data.dt);
elseif r.isSet()
    
    r.make(0,216.5e-3,2.5,r.data.dt(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, dt = %.3f ms\n',r.c.now,r.c.total,...
        r.data.dt(r.c(1))*1e3);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1);
    img = Abs_Analysis('last');
    if ~img.raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    end
    %
    % Store raw data
    %
%     r.data.img{i1,1} = img;
    r.data.w(i1,:) = img.clouds(1).becWidth;
    r.data.files{i1,1} = img.raw.files;
    %
    % Get processed data
    %
    figure(10);clf
    plot(r.data.dt(1:i1),r.data.w(1:i1,:),'o');
end