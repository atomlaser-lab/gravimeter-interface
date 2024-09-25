function Callback_Rhys_3ROI(r)
% % % Inputs
FigNum = 5;
clear TitleStuff
TitleStuff.TotalPower = '69*2';
TitleStuff.t_0 = '35';
TitleStuff.T = '5';
TitleStuff.Tau  = '30';
TitleStuff.BraggOrder = '1';

Title = append('P_{total} = ',TitleStuff.TotalPower,' mW, ', 't_0 = ',TitleStuff.t_0,' ms, ','\tau = ',TitleStuff.Tau,' us', ', T = ',TitleStuff.T,' ms', 'Bragg Order = ',TitleStuff.BraggOrder);
TitleStuff.Title = Title;
TitleStuff.SubTitle = 'Intensity Factor = 1';

Title = 'Bragg pi pulse Stability';
TitleStuff.SubTitle = 'ROI 1 is mf = 0';


Param = 25:25:400;
% Param = [0:2:360];
% Param = linspace(2.51097,2.51104,30);

PlotFactor = 1;
ParamName = 'Run';


if r.isInit()
%     r.data.param = const.randomize(Param);
    r.data.param = Param;
    r.data.PlotParam = r.data.param*PlotFactor;
    r.data.ParamName = ParamName;
    r.c.setup('var',r.data.param);

elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.param(r.c(1))).upload;
    fprintf(1,'Run  %d/%d, %s = %.3f %s\n',r.c.now,r.c.total,r.data.ParamName,r.data.param(r.c(1)));

elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = 1;
    pause(0.1 + 0.5*rand);
%     img = Abs_Analysis_GUI('last',1);    
    img = Abs_Analysis('last',1);    

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
    r.data.C1.Nsum(i1) = img.clouds(1).fit.Nsum;
    r.data.C1.xPos(i1) = img.clouds(1).fit.pos(1);
    r.data.C1.yPos(i1) = img.clouds(1).fit.pos(2);
    r.data.C1.peakOD(i1) = img.clouds(1).peakOD;

    r.data.C2.N(i1)= img.clouds(2).fit.N;
    r.data.C2.Nsum(i1)= img.clouds(2).fit.Nsum;    
    r.data.C2.xPos(i1) = img.clouds(2).fit.pos(1);
    r.data.C2.yPos(i1) = img.clouds(2).fit.pos(2);
    r.data.C2.peakOD(i1) = img.clouds(2).peakOD;

    r.data.C3.N(i1)= img.clouds(3).fit.N;
    r.data.C3.Nsum(i1)= img.clouds(3).fit.Nsum;    
    r.data.C3.xPos(i1) = img.clouds(3).fit.pos(1);
    r.data.C3.yPos(i1) = img.clouds(3).fit.pos(2);
    r.data.C3.peakOD(i1) = img.clouds(3).peakOD;    

    if r.data.C2.N(i1) < 1e3
    	r.data.C2.xPos(i1) = NaN;
    	r.data.C2.yPos(i1) = NaN;
    end
    if r.data.C1.N(i1) < 1e3
        r.data.C1.xPos(i1) = NaN;
        r.data.C1.yPos(i1) = NaN;
    end

    % % % Plot as you go
    figure(FigNum); clf
    sgtitle({['{\bf\fontsize{14}' Title '}'],TitleStuff.SubTitle});
    subplot(3,1,1)
    scatter(r.data.PlotParam(1:i1),r.data.C1.Nsum./(r.data.C1.Nsum + r.data.C2.Nsum + r.data.C3.Nsum),'r')
    hold on
    scatter(r.data.PlotParam(1:i1),r.data.C2.Nsum./(r.data.C1.Nsum + r.data.C2.Nsum + r.data.C3.Nsum),'b')
    scatter(r.data.PlotParam(1:i1),r.data.C3.Nsum./(r.data.C1.Nsum + r.data.C2.Nsum + r.data.C3.Nsum),'k')
    xlabel(ParamName)
    ylabel('N norm')

    subplot(3,1,2)
    scatter(r.data.PlotParam(1:i1),(r.data.C1.N + r.data.C2.N  + r.data.C3.N),'r')
    hold on
    scatter(r.data.PlotParam(1:i1),(r.data.C1.Nsum + r.data.C2.Nsum  + r.data.C3.Nsum),'r')
    xlabel(ParamName)
    ylabel('N Total')

    subplot(3,1,3)
    scatter(r.data.PlotParam(1:i1),r.data.C1.peakOD,'r')
    hold on
    scatter(r.data.PlotParam(1:i1),r.data.C2.peakOD,'b')
    scatter(r.data.PlotParam(1:i1),r.data.C3.peakOD,'k')
    
    xlabel(ParamName)
    ylabel('OD')
    legend('ROI 1', 'ROI 2','ROI 3')
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