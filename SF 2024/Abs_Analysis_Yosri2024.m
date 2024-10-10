function img = Abs_Analysis_GUI(varargin)


% assignin('base', 'CycleEnded', true); %this is for the app to turn off the running!
%% Parameters to check in the workspace

roiRow = [1150,1480];
roiCol = [900,1250];
atomType = 'Rb87';
tof = 217e-3; % get the 'tof' variable from the base workspace
detuning = 0; % get the 'tof' variable from the base workspace
dispOD = [0.05,.2];
plotOpt = 1;
plotROI = 1;
useFilt = 1;
filtWidth = 60e-6;
useJointFit = 0;
img_type ='drop 2';
fittype = 'twocomp2d'; %for thermal cloud or twocomp2d for BEC

Command_window_info = 1; %set to 0 if you do not want to see abs_analysis status on the command window

%% ROI of interests
roiStep = 10;
offset_region.row = [];
offset_region.col = [];


%% Imaging parameters
freqs = get_trap_freq(2,2);
imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',detuning,...
    'pixelsize',5.5e-6,'freqs',2*pi*freqs,...
    'polarizationcorrection',1,'satOD',5);

if strcmpi(img_type,'in-trap')
    rot = 180; %180
    imgconsts.magnification = 1.0285;
    imgconsts.exposureTime = 17e-6;
elseif strcmpi(img_type,'drop 2')
    rot = 90;
    imgconsts.magnification = 1.0801;
    imgconsts.exposureTime = 60e-6;
else
    warning('image system name wrong');
end


%% Load raw data
directory = 'E:\2024_2\rotten\rotten change NS in mf_1\rotten 60ms 0-10NS';
% directory = 'E:\2024_2';
% 
if nargin == 0 || (nargin == 1 && strcmpi(varargin{1},'last')) || (nargin == 2 && strcmpi(varargin{1},'last') && isnumeric(varargin{2}))
    %
    % If no input arguments are given, or the only argument is 'last', or
    % if the arguments are 'last' and a numeric array, then load the last
    % image(s).  In the case of 2 arguments, the second argument specifies
    % the counting backwards from the last image
    %
    if nargin < 2
        idx = 1;
    else
        idx = varargin{2};
    end
    args = {'files','last','index',idx,'len',3};
else
    %
    % Otherwise, parse arguments as name/value pairs for input into
    % RawImageData
    %
    if mod(nargin,2) ~= 0
        error('Arguments must occur as name/value pairs!');
    end
    args = varargin;
end
%
% This loads the raw image sets
%
raw = BinaryImageData.loadImageSets('directory',directory,'rotation',rot,args{:});

numImages = numel(raw);
plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

img = AbsorptionImage.empty;
for nn = 1:numImages
    img(nn,1) = AbsorptionImage(BinaryImageData);
end


for jj = 1:numImages
    %
    % Copy immutable properties
    %
    img(jj).constants.copy(imgconsts);
    img(jj).raw.copy(raw(jj));
    img(jj).setClouds(size(roiRow,1));
    for nn = 1:numel(img(jj).clouds)
        %       
        imgsize = size(img(jj).raw.images);
        img(jj).clouds(nn).fitdata.set('imgsize',imgsize(1:2),'roirow',roiRow(nn,:),'roiCol',roiCol(nn,:),...
            'roiStep',roiStep,'fittype',fittype,'method','y');
    end
    %
    % Create image
    %
    if size(img(jj).raw.images,3) == 2
        img(jj).makeImage;
    elseif size(img(jj).raw.images,3) == 3
        img(jj).makeImage([1,2,3]);
    else
        error('Not sure what to do here');
    end
    if useFilt
        img(jj).butterworth2D(filtWidth);
    end
    %
    % Fit clouds
    %
    if useJointFit
        img(jj).jointFit([1,2]);
    else
        img(jj).fit;
    end

    %% Plotting
    if plotOpt
        if numImages == 1
            %
            % Plot absorption data and marginal distributions when there is only 1 image
            %
            figName = 'Abs Analysis yosri Figure';
            existingFigure = findobj('Type', 'figure', 'Name', figName);
            if isempty(existingFigure)
                figure('Name', figName);
            else
                figure(existingFigure);
                clf;
            end
            img(jj).plotAllData(dispOD,plotROI);
            h = gcf;
            ax = h.Children(end);
            %             title(ax,sprintf('%s, %.2f um',ax.Title.String,img.clouds.pos(2)*1e6 - 3.3048e3));
            title(ax,sprintf('%s, %.2e',ax.Title.String,img.clouds.N));

% This for zooming on the two or three ROI and discard the rest
%             if Abs_Analysis_parameters.ROI2 == 1 || Abs_Analysis_parameters.ROI3 == 1
%                 h.Children(4).XLim = [Left Right];
%                 h.Children(4).YLim = [Top Bottom];
%             end

            %this ito set a pixel marker to locate mag trap and dipoles
%             if (pixel_marker)
%                 %show the marker
%                 plot(ax(end),ax(end).XLim,[pixel_marker(2),pixel_marker(2)],'g--')
%                 plot(ax(end),[pixel_marker(1),pixel_marker(1)],ax(end).YLim,'g--')
%             end

        else
            %
            % Plot only the absorption data in a grid when there is more than one image
            %
            %             if jj == 1
            %                 figure(3);clf;
            %                 dimSubPlot=ceil(sqrt(numImages));
            %             end
            %
            %             figure(3);
            %             subplot(dimSubPlot,dimSubPlot,jj);
            %             img(jj).plotAbsData(dispOD,plotROI);
            figure(3);clf;
            img(jj).plotAllData(dispOD,plotROI);
            pause(0.01);
        end
    end

    %% Print summaries
    [labelStr,numStr] = img(jj).labelOneROI;
    if jj == 1
        disp(labelStr);
    end
    %     for nn = 1:numel(img(jj).clouds)
    %         [labelStr,numStr] = img(jj).labelOneROI;
    disp(numStr);


% 
%     if Abs_Analysis_parameters.ROI2 == 1
%         for nn = 1:numel(img(jj).clouds)
%             R=img.clouds(2).Nsum/(img.clouds(1).Nsum+img.clouds(2).Nsum);
%         end
%         fprintf("Ratio = %.3f%%\n", R*100);
%     end
%  
%  if Abs_Analysis_parameters.ROI2 == 0
%     if strcmpi(Abs_Analysis_parameters.fittype,'twocomp2d')
%         for nn = 1:numel(img(jj).clouds)
%             Nbec=img.clouds(1).Nsum*img.clouds.becFrac;
%         end
%         fprintf("Nbec = %.4e\n", Nbec);
%     end
 end
    if Command_window_info == 1

        if useFilt == 1
            cprintf('comments', 'Camera: %s | Fit type: %s | Drop Time: %.2f ms | Detuning: %.1f MHz | Filter ON: %.1f µm\n',img_type, fittype, tof*1000, detuning, filtWidth*1e6)
        else
            cprintf('comments', 'Camera: %s | Fit type: %s | Drop Time: %.2f ms | Detuning: %.1f MHz | No Filter\n',img_type, fittype, tof*1000, detuning)
        end

    end
end

