function Callback_Rhys(r)
% r.reset empties the r.data field. to call this function, you must use
% r.reset;r.data. So, I don't think you can store data from a previous scan
% in this script

% I need to save the data then run this script


% % % Inputs
ClearImage = 1;
FigNum = 5;
% TOF = 216.5e-3;
TOF = 20e-3;

Title = 'Raman Test: TOF  ms, detuning = 0';


Param = [-15:1:10];

ParamName = 'AOM Power';
Unit = ' (Arb.)';  %do not forget to put a space before the unit

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
    r.data.becFrac(i1,:) = img.clouds.becFrac;
    r.data.N(i1,:) = img.get('N');
    r.data.T(i1,:) = img.get('T');
    r.data.Nsum(i1,:) = img.get('Nsum');
    r.data.peakOD(i1,:) = img.get('peakOD');
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:),2);
    r.data.xwidth(i1,i2) = img.clouds.gaussWidth(1);
    r.data.ywidth(i1,i2) = img.clouds.gaussWidth(2);
    r.data.BECWidthX(i1,i2) = img.clouds.becWidth(1);
    r.data.BECWidthY(i1,i2) = img.clouds.becWidth(2);
    r.data.xPos(i1,i2) = img.clouds.fit.pos(1);
    r.data.yPos(i1,i2) = img.clouds.fit.pos(2);

    % % % Plot as you go
    figure(FigNum);
    if ClearImage == 1
        clf
    end
    sgtitle(Title)
    %(x,y,YLabel,SubPlot,XName,XUnit,color,Index)
    MakePlot(r.data.param,r.data.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'r',i1)
    MakePlot(r.data.param,r.data.Nsum,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'b',i1,['Fit';'Sum'])

    MakePlot(r.data.param,r.data.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,'b',i1)

    MakePlot(r.data.param,r.data.xwidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'r',i1)
    MakePlot(r.data.param,r.data.ywidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])

    %     MakePlot(r.data.param,r.data.becFrac,'BEC Frac',4,r.data.ParamName,r.data.ParamUnits,'r',i1)
    MakePlot(r.data.param,r.data.T(:,1)*1e6,'Temperature (uK)',4,r.data.ParamName,r.data.ParamUnits,'r',i1)
    MakePlot(r.data.param,r.data.T(:,2)*1e6,'Temperature (uK)',4,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])
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