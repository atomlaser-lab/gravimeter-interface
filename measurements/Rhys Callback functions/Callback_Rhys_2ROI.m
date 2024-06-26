function Callback_Rhys_2ROI(r)
% r.reset empties the r.data field. to call this function, you must use
% r.reset;r.data. So, I don't think you can store data from a previous scan
% in this script

% I need to save the data then run this script


% % % Inputs
ClearImage = 0;
FigNum = 7;
% TOF = 216.5e-3;
TOF = 33e-3;

Title = 'Raman';
Param = (0:5:90)*pi/180;

ParamName = '4*AOM Setting (MHz)';


if r.isInit()
    r.data.param = const.randomize(Param);
%     r.data.param = (Param);
    r.data.ParamName = ParamName;
%     r.data.ParamUnits = Unit;
    r.c.setup('var',r.data.param);

elseif r.isSet()
%     r.devices.opt.detuning = r.data.param(r.c(1));
    r.make(r.devices.opt,'params',r.data.param(r.c(1)), 'tof', TOF).upload;

    fprintf(1,'Run  %d/%d, %s = %.3f %s\n',r.c.now,r.c.total,r.data.ParamName,r.data.param(r.c(1)));

elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = 1;
    pause(0.1 + 0.5*rand);
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
    r.data.C1.xPos(i1) = img.clouds(1).fit.pos(1);
    r.data.C1.yPos(i1) = img.clouds(1).fit.pos(2);
    r.data.C1.peakOD(i1) = img.clouds(1).peakOD;

    r.data.C2.N(i1)= img.clouds(2).fit.N;
    r.data.C2.xPos(i1) = img.clouds(2).fit.pos(1);
    r.data.C2.yPos(i1) = img.clouds(2).fit.pos(2);
    r.data.C2.peakOD(i1) = img.clouds(2).peakOD;

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
    scatter(r.data.param(1:i1),r.data.C1.N./(r.data.C1.N + r.data.C2.N),'r')
    hold on
    scatter(r.data.param(1:i1),r.data.C2.N./(r.data.C1.N + r.data.C2.N),'b')
    xlabel(ParamName)
    ylabel('N norm')

    subplot(3,1,2)
    scatter(r.data.param(1:i1),(r.data.C1.N + r.data.C2.N),'r')
    xlabel(ParamName)
    ylabel('N Total')

    subplot(3,1,3)
    scatter(r.data.param(1:i1),r.data.C1.peakOD,'r')
    hold on
    scatter(r.data.param(1:i1),r.data.C2.peakOD,'b')
    xlabel(ParamName)
    ylabel('OD')
    legend('ROI 1', 'ROI 2')


    % % % Plot as you go
    figure(FigNum+1);
    if ClearImage == 1
        clf
    end
    sgtitle(Title)
    subplot(3,1,1)
    scatter(r.data.param(1:i1),r.data.C1.N./(r.data.C1.N + r.data.C2.N),'r')
    ylabel('ROI1 N norm')
    xlabel(ParamName)

    subplot(3,1,2)
    scatter(r.data.param(1:i1),r.data.C2.N./(r.data.C1.N + r.data.C2.N),'b')
    xlabel(ParamName)
    ylabel('ROI2 N norm')

    subplot(3,1,3)
    scatter(r.data.param(1:i1),r.data.C1.N)
    ylabel('N total')
    xlabel(ParamName)

    subplot(3,1,3)
    scatter(r.data.param(1:i1),r.data.C1.peakOD,'r')
    hold on
    scatter(r.data.param(1:i1),r.data.C2.peakOD,'b')
    xlabel(ParamName)
    ylabel('OD')
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