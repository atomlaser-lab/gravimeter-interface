function Callback_Rhys_Bragg(r)
% r.reset empties the r.data field. to call this function, you must use
% r.reset;r.data. So, I don't think you can store data from a previous scan
% in this script

% I need to save the data then run this script


% % % Inputs
ClearImage = 0;
StorePreviousData = 1;
FigNum = 334;
% TOF = 216.5e-3;
TOF = 35e-3;

Title = 'Bragg Pulse';
Param = sort(unique((1:100)));
ParamName = 'Run Number';
Unit = ' (N)';  %do not forget to put a space before the unit


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
    r.data.files{i1,:} = img.raw.files;    
    r.data.N(i1,:,:,:) = img.get('N');
    r.data.peakOD(i1,:,:,:) = img.get('peakOD');


    % % % Plot as you go
    figure(FigNum);
    if ClearImage == 1
        clf
    end
    sgtitle(Title)
    %(x,y,YLabel,SubPlot,XName,XUnit,color,Index)
    MakePlot(r.data.param,r.data.N(:,2),'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'r',i1)
%     MakePlot(r.data.param,r.data.N(:,1),'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'r',i1,['Untransferred';'Transferred'])
    MakePlot(r.data.param,r.data.peakOD(:,2),'OD',2,r.data.ParamName,r.data.ParamUnits,'b',i1)
%     MakePlot(r.data.param,r.data.peakOD(:,1),'OD',2,r.data.ParamName,r.data.ParamUnits,'r',i1,['Untransferred';'Transferred'])

    MakePlot(r.data.param,r.data.N(:,1)./(r.data.N(:,2) + r.data.N(:,1)),'Percentage Transferred',3,r.data.ParamName,r.data.ParamUnits,'r',i1)
   
    MakePlot(r.data.param,(r.data.N(:,2) + r.data.N(:,1)),'Total Atom Number',3,r.data.ParamName,r.data.ParamUnits,'r',i1)

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