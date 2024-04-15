function Callback_Rhys_2State_MW(r)

% % % Inputs
ClearImage = 0;
FigNum = 5;
TOF = 34e-3;

Title = 'In Trap MW Transfer';
Param = -50:2:50;
ParamName = 'df (kHz)';


if r.isInit()
    r.data.freq1 = const.f_Rb_groundHFS - 315e3 + 2*Param*1e3;

    r.data.duration = Param;
    r.c.setup('var',r.data.duration);

elseif r.isSet()
    r.devices.mku.writeList(r.data.freq1(r.c(1))/2,r.data.freq1(r.c(1))/2);
    r.make(r.devices.opt, 'tof', TOF).upload;
    fprintf(1,'Run %d/%d, df = %.0f kHz\n',r.c.now,r.c.total,...
        r.data.duration(r.c(1)));

elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.5*rand);
    img = Abs_Analysis_DualState_RT('last');
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
    %
    % Store raw data
    %
    r.data.files{i1,1} = img(1).raw.files;
    %
    % Get processed data
    %
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');
    if numel(img) > 1
        r.data.N(i1,2) = r.data.N(i1,2) - r.data.N(i1,1);
        r.data.Nsum(i1,2) = r.data.Nsum(i1,2) - r.data.Nsum(i1,1);
    end
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));

    figure(FigNum);
    if ClearImage == 1
        clf;
    end
    subplot(1,2,1)
    scatter(r.data.duration(1:i1),r.data.N(1:i1,:),'filled'); %,'o'
    plot_format(ParamName,'Number','',12);

    title(' Raman frequency using fit over OD')
    grid on
    hold on;
    ylim([0,Inf]);

    subplot(1,2,2)
    scatter(r.data.duration(1:i1),r.data.Rsum(1:i1,:),'filled'); %,'o'
    hold off;
    plot_format(ParamName,'Population','',12);
    ylim([0,1])
    grid on;

    sgtitle(Title)
end