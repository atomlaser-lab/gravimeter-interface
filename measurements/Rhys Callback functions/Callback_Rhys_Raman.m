function Callback_Rhys_Raman(r)
% r.reset empties the r.data field. to call this function, you must use
% r.reset;r.data. So, I don't think you can store data from a previous scan
% in this script

% I need to save the data then run this script


% % % Inputs
ClearImage = 0;
FigNum = 334;
% TOF = 216.5e-3;
TOF = 5e-3;

Title = 'TPD';
Param = sort(unique((0:1/16:20)));
ParamName = '2*TPD';
Unit = ' (ms)';  %do not forget to put a space before the unit


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
    elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.F1.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end

    % % % Data stored
    % F = 1 atoms
    r.data.F1.files{i1,1} = img(1).raw.files;
    r.data.F1.becFrac(i1,:) = img(1).clouds.becFrac;
    r.data.F1.N(i1,:) = img(1).get('N');
    r.data.F1.T(i1,:) = img(1).get('T');
    r.data.F1.Nsum(i1,:) = img(1).get('Nsum');
    r.data.F1.peakOD(i1,:) = img(1).get('peakOD');
    r.data.F1.R(i1,:) = r.data.F1.N(i1,:)./sum(r.data.F1.N(i1,:),2);
    r.data.F1.xwidth(i1,i2) = img(1).clouds.gaussWidth(1);
    r.data.F1.ywidth(i1,i2) = img(1).clouds.gaussWidth(2);
    r.data.F1.BECWidthX(i1,i2) = img(1).clouds.becWidth(1);
    r.data.F1.BECWidthY(i1,i2) = img(1).clouds.becWidth(2);
    r.data.F1.xPos(i1,i2) = img(1).clouds.fit.pos(1);
    r.data.F1.yPos(i1,i2) = img(1).clouds.fit.pos(2);
    % F = 2 atoms
    r.data.F2.files{i1,1} = img(2).raw.files;
    r.data.F2.becFrac(i1,:) = img(2).clouds.becFrac;
    r.data.F2.N(i1,:) = img(2).get('N');
    r.data.F2.T(i1,:) = img(2).get('T');
    r.data.F2.Nsum(i1,:) = img(2).get('Nsum');
    r.data.F2.peakOD(i1,:) = img(2).get('peakOD');
    r.data.F2.R(i1,:) = r.data.F2.N(i1,:)./sum(r.data.F2.N(i1,:),2);
    r.data.F2.xwidth(i1,i2) = img(2).clouds.gaussWidth(1);
    r.data.F2.ywidth(i1,i2) = img(2).clouds.gaussWidth(2);
    r.data.F2.BECWidthX(i1,i2) = img(2).clouds.becWidth(1);
    r.data.F2.BECWidthY(i1,i2) = img(2).clouds.becWidth(2);
    r.data.F2.xPos(i1,i2) = img(2).clouds.fit.pos(1);
    r.data.F2.yPos(i1,i2) = img(2).clouds.fit.pos(2);


    % % % Plot as you go
    figure(FigNum);
    if ClearImage == 1
        clf
    end
    sgtitle(Title)
    subplot(2,1,1)
    scatter(r.data.param(1:i1),r.data.F1.N(1:i1,:)./(r.data.F1.N(1:i1,:) + r.data.F2.N(1:i1,:)),'b')
    hold on
    scatter(r.data.param(1:i1),r.data.F2.N(1:i1,:)./(r.data.F1.N(1:i1,:) + r.data.F2.N(1:i1,:)),'r')
    legend('F1','F2')
    subplot(2,1,2)
    scatter(r.data.param(1:i1),r.data.F1.N(1:i1,:))


    %     %(x,y,YLabel,SubPlot,XName,XUnit,color,Index)
    %     MakePlot(r.data.param,r.data.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'r',i1)
    %     MakePlot(r.data.param,r.data.Nsum,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'b',i1,['Fit';'Sum'])
    %
    %     MakePlot(r.data.param,r.data.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,'b',i1)
    %
    %     MakePlot(r.data.param,r.data.xwidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'r',i1)
    %     MakePlot(r.data.param,r.data.ywidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])
    %
    %     %     MakePlot(r.data.param,r.data.becFrac,'BEC Frac',4,r.data.ParamName,r.data.ParamUnits,'r',i1)
    %     %     MakePlot(r.data.param,r.data.T(1)*1e6,'Temperature (uK)',4,r.data.ParamName,r.data.ParamUnits,'b',i1)
    %     %     MakePlot(r.data.param,r.data.T(2)*1e6,'Temperature (uK)',4,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])
    %     MakePlot(r.data.param,r.data.xPos*1.0285/5.5e-6,'Pixel Number',4,r.data.ParamName,r.data.ParamUnits,'r',i1)
    %     MakePlot(r.data.param,r.data.yPos*1.0285/5.5e-6,'Pixel Number',4,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])

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