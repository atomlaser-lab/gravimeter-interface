function Callback_Rhys_MW(r)
% r.reset empties the r.data field. to call this function, you must use
% r.reset;r.data. So, I don't think you can store data from a previous scan
% in this script

% I need to save the data then run this script


% % % Inputs
ClearImage = 0;
FigNum = 5;
TOF = 36e-3;

Title = 'Microwave Transfer';
Param = sort(unique((16:0.2:20)));
% T1 = (-1.5:0.1:-3);
% T2 = (-5:0.1:-3);
% T1 = [-9.8 -9.5 -9.2];
% T2 = [-6.8, -7, -7.2];
% T3 = [-4.2 -4.5 -4.7];
% T4 = [-1.8 -2 -2.2];
% T5 = [-0.2 -0.5 -0.7];
% T6 = [0.8 1 1.2];
% T7 = [2.3 2.5 2.7];
% T8 = [5.8 5.9 5.7];
% Param = sort(unique([T1 T2 T3 T4 T5 T6 T7 T8]));
% Param = sort(unique([T1 T2]));

ParamName = 'MW Frequency';
Unit = ' (kHz)';  %do not forget to put a space before the unit


if r.isInit()
    r.data.df = Param;
    r.data.freq1 = const.f_Rb_groundHFS - 315e3 + Param*1e3;
    r.data.freq2 = const.f_Rb_groundHFS*ones(size(Param));
    
    r.data.ParamName = ParamName;
    r.data.ParamUnits = Unit;
    r.c.setup('var',r.data.df);

elseif r.isSet()

    r.devices.mku.writeList(r.data.freq1(r.c(1))/2,r.data.freq2(r.c(1))/2);
    r.make(r.devices.opt, 'tof', TOF).upload;

    fprintf(1,'Run  %d/%d, %s = %.3f %s\n',r.c.now,r.c.total,r.data.ParamName,r.data.df(r.c(1)),r.data.ParamUnits);

elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = 1;
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis_GUI('last',1);
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(img.raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end

    % % % Data stored
    r.data.files{i1,1} = img.raw.files;
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');


    r.data.C1.N(i1) = img.clouds(1).fit.N;
    r.data.C1.xPos(i1) = img.clouds(1).fit.pos(1);
    r.data.C1.yPos(i1) = img.clouds(1).fit.pos(2);

    r.data.C2.N(i1)= img.clouds(2).fit.N;
    r.data.C2.xPos(i1) = img.clouds(2).fit.pos(1);
    r.data.C2.yPos(i1) = img.clouds(2).fit.pos(2);

    if r.data.C2.N(i1) < 1e3
    	r.data.C2.xPos(i1) = NaN;
    	r.data.C2.yPos(i1) = NaN;
    end
    if r.data.C1.N(i1) < 1e3
        r.data.C1.xPos(i1) = NaN;
        r.data.C1.yPos(i1) = NaN;
    end


    % % % Plot as you go
    figure(FigNum);
    if ClearImage == 1
        clf
    end
    sgtitle(Title)
    subplot(3,1,1)
    scatter(r.data.df(1:i1),r.data.C1.N./(r.data.C1.N + r.data.C2.N),'r')
    hold on
    scatter(r.data.df(1:i1),r.data.C2.N./(r.data.C1.N + r.data.C2.N),'b')
    xlabel('df')
    ylabel('N norm')
    legend('ROI 1', 'ROI 2')

    subplot(3,1,2)
    scatter(r.data.df(1:i1),r.data.C1.N + r.data.C2.N)
    ylabel('N total')
    subplot(3,1,3)
    scatter(r.data.df(1:i1),r.data.C1.yPos,'r')
    hold on
    scatter(r.data.df(1:i1),r.data.C2.yPos,'b')
    xlabel('df')
    ylabel('Pos')
    legend('ROI 1', 'ROI 2')

end



end


function MakePlot(x,z,Ylabel,SubPlot,ParameterName,ParameterUnits,colour,Index,Legend)
subplot(2,2,SubPlot);hold on
scatter(x(1:Index),z,'color',colour)
enhformat([ParameterName,ParameterUnits], Ylabel,'small')
if exist("Legend") == 1
    legend(Legend)
end
end