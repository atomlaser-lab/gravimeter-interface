function Callback_MeasurePolarizationCorrection(r)

if r.isInit()
    r.data.image_amp = const.randomize(repmat([0.025,0.05,0.1:0.1:1],4,1));
    r.data.pulse_duration = 5e-6./r.data.image_amp;
    r.c.setup('var',r.data.image_amp);
elseif r.isSet()
    r.make(r.devices.opt,'params',[r.data.image_amp(r.c(1)),r.data.pulse_duration(r.c(1))]).upload;
    fprintf(1,'Run %d/%d, Image amplitude = %.2f, Pulse duration = %.2f us\n',r.c.now,r.c.total,r.data.image_amp(r.c(1)),r.data.pulse_duration(r.c(1))*1e6);
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1 + 0.5*rand);

    imgconsts = AtomImageConstants('Rb87','tof',r.devices.opt.tof,'detuning',0,...
        'pixelsize',5.5e-6,'exposureTime',r.data.pulse_duration(i1),'polarizationcorrection',1,'satOD',Inf,...
        'photonsPerCount',0.64,'magnification',3.3,'freqs',2*pi*[90,20,90]);

    raw = BinaryImageData.loadImageSets('directory','D:\labview-images','rotation',-90,'files','last','index',1);
    
    
    if ~raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end

    img = AbsorptionImage(raw,imgconsts);
    img.setClouds(1);
    img.offset_region.row = [];img.offset_region.col = [];
    img.clouds.fitdata.set('roiRow',1200 + 300*[-1,1],'roiCol',920 + 300*[-1,1],'roiStep',4*[1,1],'fittype','gauss2d');
    img.makeImage([1,2,3]);
    img.ODcorr = img.ODraw;
    img.fit;
    figure(10);clf;
    img.plotAllData([0,1],0);

    r.data.files{i1,1} = img.raw.files;
    r.data.N(i1,1) = img.get('N');
    r.data.OD(i1,1) = img.get('peakOD');
    I = (raw.images(:,:,2) - raw.images(:,:,3))./img.constants.satN;
    r.data.sat_param(i1,1) = mean(mean(I(img.clouds.fitdata.roiRow(1):img.clouds.fitdata.roiRow(2),img.clouds.fitdata.roiCol(1):img.clouds.fitdata.roiCol(2))));



    figure(123);clf;
    plot(r.data.sat_param,r.data.OD,'o');
    plot_format('Saturation Parameter','OD','',10);
    grid on;
    ylim([0,Inf]);

   
end


end