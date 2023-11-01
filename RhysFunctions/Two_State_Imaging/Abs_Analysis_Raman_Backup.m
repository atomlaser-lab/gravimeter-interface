function img = Abs_Analysis_Raman_Backup(varargin)
%% Inputs
% % % Plotting Options
FigNum = 1;
dispOD = [0,1];
plotOpt = 1;
plotROI = 0;
useFilt = 1;
filtWidth = 50e-6;

% % % Region of interest options

% Fit type
fittype = 'gauss2d'; %gauss2d '2comp1d'

% Choose appropriate camera properties
img_type = 'drop 1';

% ROI 1
roiStep(1) = 10;
TopRow = 0;
BottomRow = 2048;
LeftColumn = 0;
RightColumn = 1024;

% ROI 2
roiStep(2) = 10;
TopRow2 = 0;
BottomRow2 = 2048;
LeftColumn2 = 0;
RightColumn2 = 1024;

%% Automated Inputs
% F = 2 atom imaging parameters
tof = evalin('base', 'opt.tof');
detuning = evalin('base', 'opt.detuning');

% F = 1 atom imaging parameters
tof2 = tof + evalin('base', 'opt.misc.tof2');
detuning2 = evalin('base', 'opt.misc.detuning2');

%% Constant Inputs
atomType = 'Rb87';
directory = 'E:\labview-images';
DataType = 'mono8';



%% Set ROI
roiCol{1} = [LeftColumn RightColumn; LeftColumn RightColumn];
roiRow{1} = [TopRow BottomRow; TopRow BottomRow];

roiCol{2} = [LeftColumn2 RightColumn2; LeftColumn2 RightColumn2];
roiRow{2} = [TopRow2 BottomRow2; TopRow2 BottomRow2];
%% Imaging parameters
freqs = get_trap_freq(2,2);
imgconsts{1} = AtomImageConstants(atomType,'tof',tof,'detuning',detuning,...
    'pixelsize',5.5e-6,'magnification',1.0285,...
    'freqs',2*pi*freqs,'exposureTime',5e-6,...
    'polarizationcorrection',1,'satOD',5);
imgconsts{2} = AtomImageConstants(atomType,'tof',tof2,'detuning',detuning2,...
    'pixelsize',5.5e-6,'magnification',1.0285,...
    'freqs',2*pi*freqs,'exposureTime',5e-6,...
    'polarizationcorrection',1,'satOD',5);

if strcmpi(img_type,'drop 1')
    rot = 180;
    imgconsts{1}.magnification = 1.0285;
    imgconsts{1}.exposureTime = 5e-6;
    imgconsts{2}.magnification = 1.0285;
    imgconsts{2}.exposureTime = 5e-6;
elseif strcmpi(img_type,'drop 2')
    rot = 90;
    imgconsts{1}.magnification = 1.0801;
    imgconsts{1}.exposureTime = 50e-6;
    imgconsts{2}.magnification = 1.0801;
    imgconsts{2}.exposureTime = 50e-6;
end


%% Load raw data
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
raw = BinaryImageData.loadImageSets('directory',directory,'rotation',rot,args{:},'datatype',DataType);

numImages = numel(raw);

plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

img = AbsorptionImage.empty;

if size(raw(1).images,3) <= 3
    for nn = 1:numImages
        img(nn,1) = AbsorptionImage(BinaryImageData);
    end
else
    for nn = 1:numImages
        img(nn,1) = AbsorptionImage(BinaryImageData);
        img(nn,2) = AbsorptionImage(BinaryImageData);
    end
end

for kk = 1:2
    for jj = 1:numImages

        %
        % Copy immutable properties
        %
        img(jj,kk).constants.copy(imgconsts{kk});
        img(jj,kk).raw.copy(raw(jj));
        if kk == 1
            img(jj,kk).raw.images = img(jj,kk).raw.images(:,:,[1,3,4]); % atoms in F = 2 manifold
        else
            img(jj,kk).raw.images = img(jj,kk).raw.images(:,:,[2,3,4]); % atoms in F = 1 manifold
        end
        img(jj,kk).setClouds(size(roiRow{kk},1));
        for nn = 1:numel(img(jj,kk).clouds)
            img(jj,kk).clouds(nn).fitdata.set('roirow',roiRow{kk}(nn,:),'roiCol',roiCol{kk}(nn,:),...
                'roiStep',roiStep(kk),'fittype',fittype,'method','x');
        end

        %
        % Create image (i.e. construct OD)
        %
        img(jj,kk).makeImage([1,2,3]);

        if useFilt
            img(jj).butterworth2D(filtWidth);
        end
        %
        % Fit clouds
        %
        img(jj,kk).fit;

        %% Plotting original
        if plotOpt
            if numImages == 1
                %
                % Plot absorption data and marginal distributions when there is only 1 image
                %
                figure(FigNum + (kk - 1));clf;
                img(jj,kk).plotAllData(dispOD,plotROI);
                h = gcf;
                ax = h.Children(end);

                if kk == 1
                    title(ax,sprintf('%s',[ax.Title.String, ', F =2']));
                else
                    title(ax,sprintf('%s',[ax.Title.String, ', F =1']));
                end
            else
                figure(FigNum);clf;
                img(jj,l).plotAllData(dispOD,plotROI);
                pause(0.01);
            end
        end

        %% Print summaries
        [labelStr,numStr] = img(jj,kk).labelOneROI;
        if jj == 1
            disp(labelStr);
        end
        disp(numStr);

    end
end

end

