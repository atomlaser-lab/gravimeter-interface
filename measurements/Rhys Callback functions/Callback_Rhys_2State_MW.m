function Callback_Rhys_2State_MW(r)

% % % Inputs
ClearImage = 0;
FigNum = 5;
TOF = 34e-3;

Title = 'In Trap MW Transfer';
Param = -5.2:0.05:-4.8;
ParamName = 'df (kHz)';


if r.isInit()
    r.data.freq1 = const.f_Rb_groundHFS - 315e3 + 2*Param*1e3;

    r.data.Param = Param;
    r.c.setup('var',r.data.Param);

elseif r.isSet()
    r.devices.mku.writeList(r.data.freq1(r.c(1))/2,r.data.freq1(r.c(1))/2);
    r.make(r.devices.opt, 'tof', TOF).upload;
    fprintf(1,'Run %d/%d, df = %.0f kHz\n',r.c.now,r.c.total,...
        r.data.Param(r.c(1)));

elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.5*rand);
    img = Abs_Analysis_DualState_RT('last');
    if ~img(1).raw.status.ok()
        
%         Checks for an error in loading the files (caused by a missed
%         image) and reruns the last sequence
        
        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end
    
%     Store raw data
    
    r.data.files{i1,1} = img(1).raw.files;
    
%     Get processed data
    
    r.data.F2.N(i1,:) = img(1).get('N');
    r.data.F2.Nsum(i1,:) = img(1).get('Nsum');


    r.data.F1.N(i1,:) = img(2).get('N');
    r.data.F1.Nsum(i1,:) = img(2).get('Nsum');

    
    r.data.R(i1,:) = r.data.F1.N(i1,:)./(r.data.F1.N(i1,:) + r.data.F2.N(i1,:));
    r.data.Rsum(i1,:) = r.data.F1.Nsum(i1,:)./(r.data.F1.Nsum(i1,:) + r.data.F2.Nsum(i1,:));
    
    figure(FigNum);clf;
    subplot(1,2,1)
    scatter(r.data.Param(1:i1),r.data.F2.N(1:i1,:),'filled'); %,'o'
    hold on
    scatter(r.data.Param(1:i1),r.data.F1.N(1:i1,:),'filled'); %,'o'
    plot_format(ParamName,'Number','',12);
    legend('F2','F1')
    grid on
    hold on;
    ylim([0,Inf]);
    
    subplot(1,2,2)
    scatter(r.data.Param(1:i1),r.data.R(1:i1,:),'filled','color','b'); %,'o'
    hold on
    scatter(r.data.Param(1:i1),r.data.Rsum(1:i1,:),'filled','color','r'); %,'o'    
    plot_format(ParamName,'F1 Population','',12);
    ylim([0,1])
    grid on;

    sgtitle(Title)
end