function img = Abs_Analysis_DualState(varargin)
atomType = 'Rb87';
imaging_system = 'high res';
tof = 18.1e-3;
detuning = 0;
dispOD = [0,1];
plotOpt = 1;
plotROI = 0;
useFilt = 0;
filtWidth = 5e-6;
%% Set imaging region-of-interest (ROI)
%fit types
%none
% sum
% gauss1d
% twocomp1d
% tf1d
% gauss2d
% tf2d
% twocomp2d
% pixel_marker = [1251,379];

pixel_marker = [685,445];
% %for full frame
% roiRow{1} = [10,2048];
% roiCol{1} = [10,2048];
%for alternative frame
% roiRow{1} = 750 + 250*[-1,1];
% roiCol{1} = 920 + 300*[-1,1];
roiRow{1} = 840 + 300*[-1,1];
roiCol{1} = 950 + 400*[-1,1];
roiStep{1} = 3*[1,1];
fittype{1} = 'gauss2d';

roiRow{2} = 200 + 200*[-1,1];
roiCol{2} = 950 + 300*[-1,1];
roiStep{2} = 3*[1,1];
fittype{2} = 'gauss2d';
offset_region.row = [1,100];
offset_region.col = [1,100];
%% Imaging parameters
imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',detuning,...
    'pixelsize',5.5e-6,'exposureTime',40e-6,'polarizationcorrection',1,'satOD',11);
imgconsts.freqs = 2*pi*get_trap_freq(0.725,1.34);

if strcmpi(imaging_system,'low res')
    imgconsts.magnification = 0.64;
    imgconsts.photonsPerCount = 0.2644;
    image_rotation = 90;
elseif strcmpi(imaging_system,'high res')
    imgconsts.magnification = 3.5;
    imgconsts.photonsPerCount = 0.4747;
    image_rotation = -90; %-90
elseif strcmpi(imaging_system,'vertical')
%      imgconsts.magnification = 3.8; % 10ms drop
%         imgconsts.magnification = 6.16; % 16ms drop
    imgconsts.magnification = 7.8; % 37ms drop
    imgconsts.photonsPerCount = 0.4747;
    image_rotation = -90;
elseif strcmpi(imaging_system,'MOT')
    imgconsts.magnification = 0.75;
    imgconsts.photonsPerCount = 0.27;
    image_rotation = 90;
else
    warning('image system name wrong');
end
% directory = 'D:\raw-images';
directory = 'D:\labview-images';
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
    args = {'files','last','index',idx};
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
raw = BinaryImageData.loadImageSets('directory',directory,'rotation',image_rotation,args{:});
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
        img(jj,kk).constants.copy(imgconsts);
        img(jj,kk).raw.copy(raw(jj));
        if kk == 1
            img(jj,kk).raw.images = img(jj,kk).raw.images(:,:,[2,3,4]);
        else
            img(jj,kk).raw.images = img(jj,kk).raw.images(:,:,[1,3,4]);
        end
        img(jj,kk).setClouds(size(roiRow{kk},1));
        img(jj,kk).offset_region = offset_region;
        for nn = 1:numel(img(jj,kk).clouds)
            img(jj,kk).clouds(nn).fitdata.set('roirow',roiRow{kk}(nn,:),'roiCol',roiCol{kk}(nn,:),...
                'roiStep',roiStep{kk},'fittype',fittype{kk},'method','x');
        end
        %
        % Create image
        %
        img(jj,kk).makeImage([1,2,3]);
        if useFilt
            img(jj).butterworth2D(filtWidth);
        end
        %
        % Fit clouds
        %
        img(jj,kk).fit;
        %% Plotting
        if plotOpt
            if numImages == 1
                %
                % Plot absorption data and marginal distributions when there is only 1 image
                %
                figure(kk);clf;
                img(jj,kk).plotAllData(dispOD,plotROI);
                axs = get(gcf,'children');
                axs(end).Title.String = [axs(end).Title.String,', ',sprintf('F = %d, N = %.2e',kk,img(jj,kk).clouds(1).N)];
                if (pixel_marker)
                    %show the marker
                    plot(axs(end),axs(end).XLim,[pixel_marker(2),pixel_marker(2)],'g--')
                    plot(axs(end),[pixel_marker(1),pixel_marker(1)],axs(end).YLim,'g--')
                end

            else
                figure(kk);clf;
                img(jj,kk).plotAllData(dispOD,plotROI);
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