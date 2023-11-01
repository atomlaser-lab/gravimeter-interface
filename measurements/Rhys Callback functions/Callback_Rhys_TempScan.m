function Callback_Rhys_TempScan(r)
% % %  I cannot get the final plot to work (once scan is done)

% % % Plotting while you go is hard to understand and doesn't appear to
% work now. I will try and combine multiple 2D plots for visualisation

% I need to sort the data for the surface plot


% % % Inputs
FigNum = 333;
Title = 'MOT';
Param = sort(unique([(0:1:1)]));
Param1Name = '3D Coil Current';
Unit = ' (A)';  %do not forget to put a space before the unit
tof = (20:5:30)*1e-3;



if r.isInit()
    r.data.param = const.randomize(Param);
    r.data.ParamName = Param1Name;
    r.data.ParamUnits = Unit;

    r.data.tof = tof;
    r.data.ParamName2 = 'Time of Flight';
    r.data.ParamUnits2 = ' (ms)';

    r.c.setup('var',r.data.tof,r.data.param);

elseif r.isSet()
    r.make(r.devices.opt, 'tof', r.data.tof(r.c(1)),'params',r.data.param(r.c(2))).upload;

    fprintf(1,'Run  %d/%d, %s = %.3f %s, tof = %.3f \n',r.c.now,r.c.total,r.data.ParamName,r.data.param(r.c(1)),r.data.ParamUnits,r.data.tof(r.c(2)));


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

    % % % % plot
    figure(FigNum);
    sgtitle(Title)

    MakePlot(r.data.param,r.data.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'r',i1)
    MakePlot(r.data.param,r.data.Nsum,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'b',i1,['Fit';'Sum'])

    MakePlot(r.data.param,r.data.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,'b',i1)

    MakePlot(r.data.param,r.data.xwidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'r',i1)
    MakePlot(r.data.param,r.data.ywidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])

    MakePlot(r.data.param,r.data.xPos,'x Position ',4,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])
    MakePlot(r.data.param,r.data.yPos,'y Position ',4,r.data.ParamName,r.data.ParamUnits,'b',i1,['x';'y'])

    if r.c.done(1)
        % Prepare Data
        SigXSquared = r.data.xwidth.^2;
        SigYSquared = r.data.ywidth.^2;
        TOFSquared = tof.^2;

        % Fit X
        ft = fittype( 'poly1' );

        [xData, yData] = prepareCurveData(tofSquared, SigXSquared);
        [XFit, gof] = fit( xData, yData, ft );
        % Fit Y
        [xData, yData] = prepareCurveData(tofSquared, SigYSquared);
        [YFit, gof] = fit( xData, yData, ft );

        Tx(i2) = XFit.p1*const.mRb/const.kb*1e6;
        Ty(i2) = YFit.p1*const.mRb/const.kb*1e6;

        % Plot Fit
        figure(FigNum);
        subplot(2,2,3);hold on
        plot(tof,sqrt(XFit(tof.^2)),'r')
        plot(tof,sqrt(YFit(tof.^2)),'b')
    end


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