%Scanning Function
function ScanSinglePlotNumber(r)
if r.isInit()
    %Initialize run
    %% Analysis Requirements (enter this for your run)
    r.data.Detuning = 0;
%     r.data.Tdrop = 13e-3;
    
    %% Draw a box around your cloud below in second function
    r.data.TopEdge = 200;
    r.data.BottomEdge = 900;
    r.data.LeftEdge = 75;
    r.data.RightEdge = 700;
    
    % Note that to change the fit type, you must go to line 207, but then
    % your Nsum will be questionable (because it uses the offset obtained by a 2DGauss when calculating Nsum) 
    
    %% Simple Scan
     r.data.tof = (1:3:4); %this is varagin{1}
     r.data.param = 1;  %note that this varagin{2}

     %% expectged peaks Scan
%     RS_offset=-160e3;
%     peak1 = -0.7e5;
%     peak2 = -1.8e5;
%     peak3 = 0.5e5;
%     widths = 50e3;
    
    
	%r.data.param =(-27e5:0.1e5:0e5);
%     [s,d,f] = guassian_linespace(200,-250e3,150e3,[peak1 widths 1],[peak2 widths 1],[peak3 widths 1]);
%     figure(99);
%     subplot(1,2,1);
%     plot(f,d);
%     subplot(1,2,2);
%     plot(s);
%     r.data.param2 = ones(1,50)*49653.6;
         
    %% shuffel!
%     r.data.param =r.data.param(randperm(length(r.data.param)));
    
    %% Counter and Callback
    r.data.count = RolloverCounter([numel(r.data.tof),numel(r.data.param)]);
    r.numRuns = r.data.count.total;
    
    r.makerCallback = @MOT;
    
elseif r.isSet()
    %Increment counter
    if r.currentRun ~= 1
        r.data.count.increment;
    end
    
    %Build/upload/run sequence
    r.make(r.data.tof(r.data.count.idx(1)),r.data.param(r.data.count.idx(2)));
    r.upload;
    %Print information about current run
    fprintf(1,'Run %d/%d, MOT, param2 : %.3f \n',r.currentRun,r.numRuns,r.data.param2(r.data.count.idx(2)));

elseif r.isAnalyze()
    % Make shorthand variables for indexing
    nn = r.currentRun;
    i1 = r.data.count.idx(1);
    i2 = r.data.count.idx(2);
    pause(1.0); %Wait for other image analysis program to finish with files
%     i2=1;
    %Analyze image data from last image
    [c,Nsum] = Abs_Analysis_Internal(r.data.tof(nn),r.data.Detuning,r.data.TopEdge,r.data.BottomEdge,r.data.LeftEdge,r.data.RightEdge); %see top
    r.data.files{i1,i2} = {c.raw.files(1).name,c.raw.files(2).name};

    %% Plottable Data
    r.data.N(i1,i2) = c.N;
    r.data.Nsum(i1,i2) = Nsum;
%     r.data.Nth(i1,i2) = c.N.*(1-c.becFrac);
%     r.data.Nbec(i1,i2) = c.N.*c.becFrac;
%     r.data.F(i1,i2) = c.becFrac;
    r.data.xw(i1,i2) = c.gaussWidth(1);
%     r.data.x0(i1,i2) = c.pos(1);
    r.data.yw(i1,i2) = c.gaussWidth(2);
%     r.data.y0(i1,i2) = c.pos(2);
%     r.data.OD(i1,i2) = c.fitdata.params.gaussAmp(1);
    r.data.OD(i1,i2) = c.peakOD;
%     r.data.Tx(i1,i2) = c.T(1);
%     r.data.Ty(i1,i2) = c.T(2);
    
    %% Catch bad fits
    %need to make all variables nil
    if c.N > 1e10
       r.data.N(i1,i2) = NaN;
    end    
    
    %% Plot 2 variables on same graph(eg hold variables constant and measure fluctuation)
    data_y = r.data.N;
    data_x = (1:1:length(r.data.N));
    figure(23)
    plot(data_x,data_y,'o-');
    hold on;
    plot(data_x,r.data.Nsum(1:i1,:),'sq-');
    hold off;
    title(char(datetime))

    
    %% Plot (y1,x) on subplot 1 and (y2,x) on subplot 2
    
%    [data_x, sortIdx] = sort(r.data.param2(1:i2));
%    data_y = r.data.N;
%    data_y = data_y(sortIdx);
%    data_y_2 = r.data.OD;
%    data_y_2 = data_y_2(sortIdx);
%    
%    figure(24);%clf;
%    subplot(1,2,1);
%    cla;
%    plot(data_y,'X-');
%    subplot(1,2,2);
%    plot(data_y_2,'X-');
%    title(char(datetime))
   
    %% Plot (y,param2) for a given param1( on subplot 1 and (y,param2) on subplot 2 where each line is a different param1
    
%    [data_x, sortIdx] = sort(r.data.param(1:i1));
%    data_y = r.data.N(:,i2);
%    data_y = data_y(sortIdx);
% 
%     figure(24);
%     subplot(1,2,1);
%     cla;
%     plot(data_x,data_y,'o-');
%     subplot(1,2,2);
%     if length(data_x) == length(r.data.param)
%         if i2 == 1
%             cla;
%         end
%         s = {};
%         for jj = 1:i2
%             s{jj} = sprintf('%.3f',r.data.param2(jj));
%         end
%         plot(data_x,data_y,'o-')
%         hold on
%         legend(s);
%     end
%     title(char(datetime))
    
    %% Surface Plot of (y,param1,param2)
    
%     clf(figure(23))
%     if r.currentRun == r.numRuns
%         figure(23);
%         [x,y] = meshgrid(r.data.param,r.data.param2);
%         surf(x,y,r.data.N);
%     end
%    title(char(datetime))

end %end if statement

end %end function














%% abs analysis function called above:
function [cloud,Nsum] = Abs_Analysis_Internal(tof,detuning,TopEdge,BottomEdge,LeftEdge,RightEdge)
atomType = 'Rb87';

col1 = 'b.-';
col2 = 'r--';
dispOD = [0,0.5];
plotROI = 0;

numPixels = abs((TopEdge-BottomEdge)*(LeftEdge-RightEdge));
if numPixels > 200^2 && numPixels < 400^2
    stepSize = 2;
elseif numPixels >= 400^2 && numPixels < 600^2
    stepSize = 5;
elseif numPixels >= 600^2
    stepSize = 10;
else
    stepSize = 1;
end

fitdata = AtomCloudFit('roiRow',[TopEdge,BottomEdge],...%was 201
                       'roiCol',[LeftEdge,RightEdge],...
                       'roiStep',stepSize,...
                       'fittype','gauss2d');    %Options: none, gauss1d, twocomp1d, bec1d, gauss2d, twocomp2d, bec2d, sum

imgconsts = AtomImageConstants(atomType,'exposureTime',100e-6,...
            'pixelsize',5.5e-6,'magnification',0.25,...
            'freqs',2*pi*[40,23,8],'detuning',detuning,... %set detuning here
            'polarizationcorrection',1.5,'satOD',3);

directory = 'D:\RawImages\2020\12December\';

%% Load raw data
raw = RawImageData('filenames','last','directory',directory);

%% Analyze data
cloud = AbsorptionImage(raw,imgconsts,fitdata);
cloud.makeImage;
cloud.fit([],tof,'y');  %'y' indicates the marginal distribution to use for calculating number of atoms.  Can also be 'x' or 'xy'

offset = fitdata.params.offset;               %Determine offset on image and uses this when summing the pixels for atom number. This reduces the effects of linear ofsets
Nsum = cloud.sum(offset);
% 
% fitdata.image = fitdata.image - offset;     %Subtracts offset from image
% fitdata.fittype = 'sum';                    %Change fit type to 'sum'
% cloud.fit([],tof,'y');                      %Re-'fit' the data
% Nsum = cloud.N;

%% Plot image and fits
figure(3);clf;
cloud.plotAllData(dispOD,col1,col2,plotROI);

%% Print labels
[labelStr,numStr] = cloud.labelOneROI;
disp(labelStr);
disp(numStr);

end