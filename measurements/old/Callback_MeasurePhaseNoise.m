function Callback_MeasurePhaseNoise(r)

if r.isInit()
    r.data.run = 1:50;
    r.data.T = [25 5 15 1 10]*1e-3;
%     r.data.phase = mod([16.67,38.5,148.78,-12.815,-16.45]*2,360);
r.data.phase = [97.348,39.027,92.368,44.509,45.112];

    r.c.setup('var',r.data.run,r.data.T);
elseif r.isSet()
    
    r.make(0,216.6e-3,1.175,0.235,r.data.phase(r.c(2)),0,r.data.T(r.c(2)));
%     r.make(0,216.6e-3,1.425,0.241,r.data.phase(r.c(1)),0,r.data.Tint(r.c(2)))
    r.upload;
    fprintf(1,'Run %d/%d, T = %.2f ms\n',r.c.now,r.c.total,...
        r.data.T(r.c(2))*1e3);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.1);
    img = Abs_Analysis('last',1);
    if ~img.raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return
    end
    %
    % Store raw data
    %
    r.data.files{i1,i2} = {img.raw.files(1).name,img.raw.files(2).name};
    %
    % Get processed data
    %
    r.data.N(i1,i2,:) = reshape(img.get('N'),[1,1,2]);
    r.data.Nsum(i1,:) = reshape(img.get('Nsum'),[1,1,2]);
    r.data.R(i1,i2) = r.data.N(i1,1)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,i2) = r.data.Nsum(i1,1)./sum(r.data.Nsum(i1,:));
    
    figure(98);
    subplot(1,2,1);
    plot(r.data.run(1:i1),r.data.R(1:i1,i2),'o-');
    hold on;
    plot(r.data.run(1:i1),r.data.Rsum(1:i1,i2),'sq-');
    hold off;
    plot_format('Run number','N_{rel}','',12);
    subplot(1,2,2);
    if r.c.done(1)
        cla;
        for mm = 1:i2
            plot(r.data.run,r.data.Rsum(:,mm),'o-');
            hold on;
            str{mm} = sprintf('Time = %d ms',round(1e3*r.data.T(mm)));
        end
        plot_format('Run number','N_{rel}','',12);
        legend(str);
    end
    pause(0.01);
    
end