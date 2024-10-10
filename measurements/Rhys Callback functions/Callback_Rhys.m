function Callback_Rhys(r)
% % % Inputs
ClearImage = 1;
FigNum = 5;

Title = 'Stability NO MW';
% Param = 30:30:300;
Param = 3:0.5:14;

ParamName = 'Run';
PlotFactor = 1;

if r.isInit()
    r.data.param = const.randomize(Param);
%     r.data.param = (Param);
    r.data.plotparam = r.data.param*PlotFactor;
    r.data.ParamName = ParamName;
    r.c.setup('var',r.data.param);

elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.param(r.c(1))).upload;

    fprintf(1,'Run  %d/%d, %s = %.3f %s\n',r.c.now,r.c.total,r.data.ParamName,r.data.param(r.c(1)));

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
    r.data.Tx(i1,:) = img.clouds.T(1);
    r.data.Ty(i1,:) = img.clouds.T(2);
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
    figure(FigNum);clf
    subplot(3,1,1)
    scatter(r.data.plotparam(1:i1),r.data.N(1:i1),'filled');
    plot_format(ParamName,'N_{F = i}','',12);
    grid on

    subplot(3,1,2)
    scatter(r.data.plotparam(1:i1),r.data.peakOD(1:i1,1),'filled');
    plot_format(ParamName,'Peak OD','',12);

    subplot(3,1,3)
    scatter(r.data.plotparam(1:i1),r.data.Tx(1:i1),100,'b');
    hold on
    scatter(r.data.plotparam(1:i1),r.data.Ty(1:i1),100,'r');
    plot_format(ParamName,'Temp','',12);
    legend('x','y')
    grid on;
    
    sgtitle(Title)
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