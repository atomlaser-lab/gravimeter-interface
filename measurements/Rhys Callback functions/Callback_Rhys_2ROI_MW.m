function Callback_Rhys_2ROI_MW(r)
FigNum = 5;
Title = 'Out of Trap MW Pulse 1';
% Param = ([-10:1:20 -40:1:-21 21:1:40 -60:1:-41 41:1:60 -80:1:-61 61:1:80]);
% Param = unique([-1.5:0.5:2]) -292.75;
Param = unique([-0.6:0.2:0.6]) -287.2;

ParamName = 'MW Freq';


if r.isInit()
    r.data.df = Param;
%     r.data.freq1 = const.f_Rb_groundHFS - 434.55*1e3 + Param*1e3; % 2 power supply, EW normal
%     r.data.freq1 = const.f_Rb_groundHFS - 287.2*1e3 + Param*1e3; %1 power supply, EW flipped, 11 ms tof

    r.data.freq1 = const.f_Rb_groundHFS + Param*1e3; 


    r.data.freq2 = const.f_Rb_groundHFS;
    r.data.ParamName = ParamName;
    r.c.setup('var',r.data.df);

elseif r.isSet()
    if numel(r.data.freq1) > 1
        r.devices.mku.writeList(r.data.freq1(r.c(1))/2,r.data.freq2/2);
    else
        r.devices.mku.writeList(r.data.freq1/2,r.data.freq2(r.c(1))/2);
    end
    r.make(r.devices.opt).upload;
    fprintf(1,'Run %d/%d, T = %.0f us\n',r.c.now,r.c.total,...
        r.data.df(r.c(1)));

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
    figure(FigNum);clf
    sgtitle(Title)
    subplot(3,1,1)
    scatter(r.data.df(1:i1),r.data.C1.N./(r.data.C1.N + r.data.C2.N),'r')
    hold on
    scatter(r.data.df(1:i1),r.data.C2.N./(r.data.C1.N + r.data.C2.N),'b')
    xlabel(ParamName)
    ylabel('N norm')

    subplot(3,1,2)
    scatter(r.data.df(1:i1),(r.data.C1.N),'r')
    hold on
    scatter(r.data.df(1:i1),(r.data.C2.N),'b')
    scatter(r.data.df(1:i1),(r.data.C1.N + r.data.C2.N),'k')
    legend('total','ROI 1','ROI 2')
    xlabel(ParamName)
    ylabel('N')

    subplot(3,1,3)
    scatter(r.data.df(1:i1),r.data.C1.peakOD,'r')
    hold on
    scatter(r.data.df(1:i1),r.data.C2.peakOD,'b')
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