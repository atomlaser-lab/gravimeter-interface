figure(333);clf
% if r.c.final(2) > 1
if size(r.c,1) > 1 || size(r.c,2) > 1
    NumScannedVariables = 2;
else
    NumScannedVariables = 1;
end

if NumScannedVariables == 1

%     MakePlot2D(r.data.param,r.data.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'r')
%     MakePlot2D(r.data.param,r.data.Nsum,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'b',['Fit';'Sum'])
% 
%     MakePlot2D(r.data.param,r.data.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,'b')
% 
%     MakePlot2D(r.data.param,r.data.xwidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'r')
%     MakePlot2D(r.data.param,r.data.ywidth,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'b',['x';'y'])
% 
%     MakePlot2D(r.data.param,r.data.T(1),'Temperature (uK)',4,r.data.ParamName,r.data.ParamUnits,'r')
%     MakePlot2D(r.data.param,r.data.T(2),'Temperature (uK)',4,r.data.ParamName,r.data.ParamUnits,'b',['x';'y'])

    MakePlot2D(r.data.param,r.data.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'r')
    MakePlot2D(r.data.param,r.data.Nsum,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,'b',['Fit';'Sum'])

    MakePlot2D(r.data.param,r.data.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,'b')

    MakePlot2D(r.data.param,r.data.BECWidthY,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'r')
    MakePlot2D(r.data.param,r.data.BECWidthY,'Width (mm)',3,r.data.ParamName,r.data.ParamUnits,'b',['x';'y'])

    MakePlot2D(r.data.param,r.data.becFrac,'Temperature (uK)',4,r.data.ParamName,r.data.ParamUnits,'r')

elseif NumScannedVariables == 2
    % % % %     Sort Data
    Sorted_r_data = SortData(r);
    % realocate
    r.data = Sorted_r_data;
    % % % %     Plot data


    MakePlot3D(r.data.param,r.data.param2,r.data.N,'Atom Number',1,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2)
    MakePlot3D(r.data.param,r.data.param2,r.data.peakOD,'OD',2,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2)
    MakePlot3D(r.data.param,r.data.param2,r.data.xwidth,'X Width (mm)',3,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2)
    MakePlot3D(r.data.param,r.data.param2,r.data.ywidth,'Y Width (mm)',4,r.data.ParamName,r.data.ParamUnits,r.data.ParamName2,r.data.ParamUnits2)
end


function MakePlot2D(x,z,Ylabel,SubPlot,ParameterName,ParameterUnits,colour,Legend)
subplot(2,2,SubPlot);hold on
scatter(x,z,'color',colour)
enhformat([ParameterName,ParameterUnits], Ylabel,'small')
if exist("Legend") == 1
    legend(Legend)
end
end

function MakePlot3D(x,y,z,zLabel,SubPlot,XName,Xunit,YName,YUnit)
subplot(2,2,SubPlot)
surf(x,y,z.')
shading interp
view(90,90)
hold on
contour(x,y,z.')

XLABEL = [XName, Xunit];
YLABEL = [YName, YUnit];
% zLabel(zLabel,'latex','FontSize',14)

xlabel(XLABEL,'Interpreter','latex','FontSize',14)
ylabel(YLABEL,'Interpreter','latex','FontSize',14)
c = colorbar;
c.Label.String = zLabel;
c.Label.Interpreter = 'latex';
c.Label.FontSize = 14;
end


function tempy = SortData(r)

% % % Keep unchanged stuff
tempy.ParamName = r.data.ParamName;
tempy.ParamUnits = r.data.ParamUnits;
tempy.ParamName2 = r.data.ParamName2;
tempy.ParamUnits2 = r.data.ParamUnits2;

% % % %     Sort Data
Sorted1 = sort(r.data.param);
Sorted2 = sort(r.data.param2);

tempy.param = Sorted1;
tempy.param2 = Sorted2;

% % % % I don't knoiw if this will not work for double ups
for ii = 1:length(r.data.param)
    if ii == length(r.data.param)
        5
    end
    for jj = 1:length(r.data.param2)
        % Position of sorted variable
        Pos1 = find(Sorted1(ii) - r.data.param == 0);
        Pos2 = find(Sorted2(jj) - r.data.param2 == 0);
        % create sorted variables
        tempy.files(ii,jj) = r.data.files(Pos1,Pos2);
        tempy.becFrac(ii,jj) = r.data.becFrac(Pos1,Pos2);
        tempy.N(ii,jj) = r.data.N(Pos1,Pos2);
        tempy.T(ii,jj,:) = r.data.T(Pos1,Pos2,:);
        tempy.Nsum(ii,jj) = r.data.Nsum(Pos1,Pos2);
        tempy.peakOD(ii,jj) = r.data.peakOD(Pos1,Pos2);
        tempy.R(ii,jj) = r.data.R(Pos1,Pos2);
        tempy.xwidth(ii,jj) = r.data.xwidth(Pos1,Pos2);
        tempy.ywidth(ii,jj) = r.data.ywidth(Pos1,Pos2);
        tempy.BECWidthX(ii,jj) = r.data.BECWidthX(Pos1,Pos2);
        tempy.BECWidthY(ii,jj) = r.data.BECWidthY(Pos1,Pos2);
        tempy.xPos(ii,jj) = r.data.xPos(Pos1,Pos2);
        tempy.yPos(ii,jj) = r.data.yPos(Pos1,Pos2);
    end
end
end