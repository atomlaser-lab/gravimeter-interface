function Callback_MeasureMWTransferStability(r)

if r.isInit()
    r.c.setup(Inf);
elseif r.isSet()
    
    r.make(r.devices.opt);
    r.upload;
    %
    % These commands are for list-mode operation
    %
    fprintf(1,'Run %d/%d\n',r.c.now,r.c.total);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1);
    img = Abs_Analysis_old('last');
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
    r.data.files{i1,1} = img.raw.files;
    %
    % Get processed data
    %
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));

%     [~,N,dout] = FMI_Analysis;
%     r.data.N(i1,:) = [N.N1,N.N2];
%     r.data.R(i1,1) = N.R;
%     r.data.d{i1,1} = dout;
    
    
    figure(98);clf;
    subplot(1,2,1)
    plot(1:i1,r.data.R(1:i1,:),'o-');
    plot_format('Run','Population','',12);
    ylim([0,1]);
%     title('')
    grid on
    hold on;
    
    subplot(1,2,2)
    plot(1:i1,r.data.N(1:i1,:),'o-');
    plot_format('Run','Number','',12);
%     h = legend('m = -1','m = 0','total');
%     set(h,'Location','West');
%     title('Raman frequency using fit over OD')
    grid on
    hold on;

end