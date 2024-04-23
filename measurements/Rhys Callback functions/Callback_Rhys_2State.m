function Callback_Rhys_2State(r)

% % Inputs
ClearImage = 1;
FigNum = 5;
TOF = 30e-3;

Title = 'Raman:P = 0.01, \tau = 5 us, t_0 = 0.5 ms';

% T1 = -20 + (-70:5:70)*1e-3;
% T2 = 20 + (0:10:100)*1e-3;
% Param = [T1 T2];
Param = 20 + (-20:2:20)*1e-3;

PlotParam = Param;
% ParamName = 'Pulse Two Phase / \pi';
ParamName = '\delta [MHz]';



if r.isInit()
    r.data.Param = Param;
    r.data.PlotParam = PlotParam;
    r.c.setup('var',r.data.Param);

elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.Param(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, T = %.0f us\n',r.c.now,r.c.total,...
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
    scatter(r.data.PlotParam(1:i1),r.data.F2.N(1:i1,:),'filled'); %,'o'
    hold on
    scatter(r.data.PlotParam(1:i1),r.data.F1.N(1:i1,:),'filled'); %,'o'
    plot_format(ParamName,'Number','',12);
    legend('F2','F1')
    grid on
    hold on;
    ylim([0,Inf]);
    
    subplot(1,2,2)
    scatter(r.data.PlotParam(1:i1),r.data.R(1:i1,:),'filled','color','b'); %,'o'
    hold on
    scatter(r.data.PlotParam(1:i1),r.data.Rsum(1:i1,:),'filled','color','r'); %,'o'
    legend('Fit','Sum')
    hold off;
    plot_format(ParamName,'F1 Population','',12);
    ylim([0,1])
    grid on;

    sgtitle(Title)
end