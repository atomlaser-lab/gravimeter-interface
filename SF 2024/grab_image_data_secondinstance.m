function ydata = grab_image_data_secondinstance(directory, varargin)
    % Check if a directory is provided, otherwise use the default directory
    if nargin < 1 || isempty(directory)
        error('A directory must be provided.');
    end
atomType = 'Rb87';
tof = 25e-3;

dispOD = [0,0.5];
plotOpt = 0;
plotROI = 0;

%% Imaging Second spot
roiRow = [1,2048];
roiCol = 1040 + 200*[-1,1];
roiStep = 1;
img_type = 'drop 2';
%% Imaging parameters
freqs = get_trap_freq(2,2);
imgconsts = AtomImageConstants(atomType,'tof',tof,'detuning',0,...
            'pixelsize',5.5e-6,'magnification',1.0285,...
            'freqs',2*pi*freqs,'exposureTime',5e-6,...
            'polarizationcorrection',1,'satOD',5);
 
if strcmpi(img_type,'drop 1')
    rot = 180;
    imgconsts.magnification = 1.0285;
    imgconsts.exposureTime = 5e-6;
elseif strcmpi(img_type,'drop 2')
    rot = 90;
    imgconsts.magnification = 1.0801;
    imgconsts.exposureTime = 50e-6;
end  

% directory = 'Z:\Gravy\2024\last\accelerometry\two-pulse interferometer consistency';

%% Load raw data
args = parse_arguments(varargin{:});
%
% This loads the raw image sets
%
raw = BinaryImageData.loadImageSets('directory',directory,'rotation',rot,'load_images',false,args{:});

numImages = numel(raw);
plotOpt = plotOpt || numImages == 1;    %This always enables plotting if only one image is analyzed

img = AbsorptionImage(BinaryImageData);

ydata = zeros(max(roiRow) - min(roiRow) + 1,numImages);

for jj = 1:numImages
    %
    % Copy immutable properties
    %
    img.constants.copy(imgconsts);
    img.raw.copy(raw(jj));
    img.raw.load('filenames',[],'directory',directory,'rotation',rot);
    img.setClouds(size(roiRow,1));
    imgsize = size(img.raw.images);
    img.clouds.fitdata.set('imgsize',imgsize(1:2),'roirow',roiRow,'roiCol',roiCol,...
        'roiStep',roiStep,'fittype','none','method','y');
    %
    
    % Create image
    %
    if size(img.raw.images,3) == 2
        img.makeImage;
    elseif size(img.raw.images,3) == 3
        img.makeImage([1,2,3]);
    else
        error('Not sure what to do here');
    end
    img.fit;
        
    %% Plotting
%     if plotOpt
%         figure(3);clf;
%         img.plotAbsData(dispOD,plotROI);
%         hold on
%         img.plotROI;
%         hold off
%     end
    
    ydata(:,jj) = img.clouds.fitdata.ydata;
    
end

ydata = ydata.';

end

function args = parse_arguments(varargin)

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

end