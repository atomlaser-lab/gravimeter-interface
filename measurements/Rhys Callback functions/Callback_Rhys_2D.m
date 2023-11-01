function Callback_Rhys_2D(r)
% % %  I cannot get the final plot to work (once scan is done)

% % % Plotting while you go is hard to understand and doesn't appear to
% work now. I will try and combine multiple 2D plots for visualisation

% I need to sort the data for the surface plot


% % % Inputs
FigNum = 333;
% TOF = 216.5e-3;
TOF = 10e-3;

Title = 'MOT';
Param = sort(unique([(0.5:0.5:1)]));
Param1Name = '3D Coil Current';
Unit = ' (A)';  %do not forget to put a space before the unit

Param2 = sort(unique([(-40:5:-0)]));
Param2Name = 'Imaging Voltage';
Unit2 = ' (V)';  %do not forget to put a space before the unit



if r.isInit()
    r.data.param = const.randomize(Param);
    r.data.ParamName = Param1Name;
    r.data.ParamUnits = Unit;

    r.data.param2 = const.randomize(Param2);
    r.data.ParamName2 = Param2Name;
    r.data.ParamUnits2 = Unit2;

    r.c.setup('var',r.data.param,r.data.param2);

elseif r.isSet()
    r.make(r.devices.opt,'params',[r.data.param(r.c(1)),r.data.param2(r.c(2))], 'tof', TOF).upload;
    %     fprintf(1,'Run  %d/%d, %s = %.3f
    %     %s\n',r.c.now,r.c.total,r.data.ParamName,r.data.param(r.c(1)),r.data.ParamUnits); % old
    fprintf(1,'Run %d/%d, Param = %.3f, Param2 = %.3f\n',r.c.now,r.c.total,...
        r.data.param(r.c(1)),r.data.param2(r.c(2)));

elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
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

    r.data.files{i1,i2} = img.raw.files;
    r.data.becFrac(i1,i2,:) = img.clouds.becFrac;
    r.data.N(i1,i2,:) = img.get('N');
    r.data.T(i1,i2,:) = img.get('T');
    r.data.Nsum(i1,i2,:) = img.get('Nsum');
    r.data.peakOD(i1,i2,:) = img.get('peakOD');
    r.data.R(i1,i2,:) = r.data.N(i1,i2,:)./sum(r.data.N(i1,i2,:),2);
    r.data.xwidth(i1,i2) = img.clouds.gaussWidth(1);
    r.data.ywidth(i1,i2) = img.clouds.gaussWidth(2);
    r.data.BECWidthX(i1,i2) = img.clouds.becWidth(1);
    r.data.BECWidthY(i1,i2) = img.clouds.becWidth(2);
    r.data.xPos(i1,i2) = img.clouds.fit.pos(1);
    r.data.yPos(i1,i2) = img.clouds.fit.pos(2);


    if r.c(1) == r.c.final(1) && r.c(2) == r.c.final(2)
        figure(FigNum);clf
        sgtitle(Title)
        %  MakePlot(x,y,z,zLabel,SubPlot,XName,Xunit,YName,YUnit,Index1,Index2)
        MakePlot(r.data.param,r.data.param2,r.data.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2,i1,i2)
        MakePlot(r.data.param,r.data.param2,r.data.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2,i1,i2)
        MakePlot(r.data.param,r.data.param2,r.data.xwidth,'X Width (mm)',3,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2,i1,i2)
        MakePlot(r.data.param,r.data.param2,r.data.ywidth,'Y Width (mm)',4,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2,i1,i2)

    else

        figure(FigNum)
        if r.c(1) == 1 && r.c(2) == 1
            clf
        end
%         PlotAsYouGo(r.data.param,r.data.param2,r.data.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2,i1,i2)
% 
%         PlotAsYouGo(r.data.param,r.data.param2,r.data.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2,i1,i2)
% 
%         PlotAsYouGo(r.data.param,r.data.param2,r.data.xwidth,'X Width (mm)',3,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2,i1,i2)
% 
%         PlotAsYouGo(r.data.param,r.data.param2,r.data.ywidth,'Y Width (mm)',4,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2,i1,i2)
    end
end

end

function MakePlot(x,y,z,zLabel,SubPlot,XName,Xunit,YName,YUnit,Index1,Index2)
subplot(2,2,SubPlot)
surf(x(1:Index1),y(1:Index2),z)
shading interp
view(90,90)
hold on
contour(x,y,z)

XLABEL = [XName, Xunit];
YLABEL = [YName, YUnit];

xlabel(XLABEL,'Interpreter','latex','FontSize',14)
ylabel(YLABEL,'Interpreter','latex','FontSize',14)
c = colorbar;
c.Label.String = zLabel;
c.Label.Interpreter = 'latex';
c.Label.FontSize = 14;
end

function PlotAsYouGo(x,y,z,zLabel,SubPlot,XName,Xunit,YName,YUnit,Index1,Index2)
subplot(2,2,SubPlot); hold on
scatter3(z,x(1:Index1),y(1:Index2))
view(45,45)
XLABEL = [XName, Xunit];
YLABEL = [YName, YUnit];

xlabel(XLABEL,'Interpreter','latex','FontSize',14)
ylabel(YLABEL,'Interpreter','latex','FontSize',14)
ylabel(zLabel,'Interpreter','latex','FontSize',14)
end