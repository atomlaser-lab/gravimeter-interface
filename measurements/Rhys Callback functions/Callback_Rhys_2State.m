function Callback_Rhys_2State(r)
% r.reset empties the r.data field. to call this function, you must use
% r.reset;r.data. So, I don't think you can store data from a previous scan
% in this script

% I need to save the data then run this script


% % % Inputs
ClearImage = 1;
FigNum = 6;
TOF = 34e-3;

Title = 'Raman: AOM Power = 1, Pulse = 50 us, TOF = 16.5 ms';

Param = sort(unique(((-0.1:0.02:0.2))));

ParamName = 'Run Number';
Unit = ' ';  %do not forget to put a space before the unit
% ParamName = 'NS Bias Amp';
% Unit = ' (V)';  %do not forget to put a space before the unit

if r.isInit()
%     r.data.param = const.randomize(Param);
    r.data.param = (Param);
    r.data.ParamName = ParamName;
    r.data.ParamUnits = Unit;
    r.c.setup('var',r.data.param);

elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.param(r.c(1)), 'tof', TOF).upload;

    fprintf(1,'Run  %d/%d, %s = %.3f %s\n',r.c.now,r.c.total,r.data.ParamName,r.data.param(r.c(1)),r.data.ParamUnits);

elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = 1;
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis_Raman('last',1);
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

    % % % Data stored
    r.data.files{i1,1} = img(1).raw.files;
    % F = 2 atom data
    r.data.C1.N(i1) = img(1).clouds.N;
    r.data.C1.peakOD(i1) = img(1).clouds.peakOD;
%     r.data.xwidth(i1,i2) = img.clouds.gaussWidth(1);
%     r.data.ywidth(i1,i2) = img.clouds.gaussWidth(2);
%     r.data.xPos(i1,i2) = img.clouds.fit.pos(1);
%     r.data.yPos(i1,i2) = img.clouds.fit.pos(2);
%     r.data.T(i1,:) = img.get('T');
    % Total atom data
    r.data.C2.N(i1) = img(2).clouds.N;
    r.data.C2.peakOD(i1) = img(2).clouds.peakOD;

    % % % Plot as you go
    figure(FigNum);
    if ClearImage == 1
        clf
    end
    sgtitle(Title)
    %(x,y,YLabel,SubPlot,XName,XUnit,color,Index)
    subplot(2,1,1); hold on
    scatter(r.data.param(1:i1),r.data.C2.N./r.data.C1.N,'color','r')
    xlabel(ParamName)
    ylabel('N_F2/N_Total')
    subplot(2,1,2); hold on
    scatter(r.data.param(1:i1),r.data.C1.peakOD,'color','r')
    scatter(r.data.param(1:i1),r.data.C2.peakOD,'color','b')
    legend('')

%     MakePlot(r.data.param,r.data.C1.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'r',i1)
%     MakePlot(r.data.param,r.data.C2.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'b',i1,['F2';'Total'])
% 
%     MakePlot(r.data.param,r.data.C1.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,'r',i1)
%     MakePlot(r.data.param,r.data.C2.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,'b',i1,['F2';'Total'])
% 
%     MakePlot(r.data.param,r.data.xwidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'r',i1)
%     MakePlot(r.data.param,r.data.ywidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])
% 
%     MakePlot(r.data.param,r.data.T(:,1)*1e6,'Temperature (uK)',4,r.data.ParamName,r.data.ParamUnits,'r',i1)
%     MakePlot(r.data.param,r.data.T(:,2)*1e6,'Temperature (uK)',4,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])




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