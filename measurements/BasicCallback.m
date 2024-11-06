function BasicCallback(r)

if r.isInit()
    r.data.param = 1:200; % Detuning
    r.c.setup('var',r.data.param);
elseif r.isSet()
%     r.make(r.devices.opt,'params',r.data.param(r.c(1))).upload;
%     r.make(r.devices.opt).upload;
    fprintf(1,'Run %d/%d, Param = %.3f\n',r.c.now,r.c.total,r.data.param(r.c(1)));
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.25*rand);
    img = Abs_Analysis('last',1);
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end
    
    r.data.files{i1,1} = img(1).raw.files;
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');
    r.data.becFrac(i1,:) = img.get('becFrac');
    r.data.OD(i1,:) = img.get('peakOD');
%     r.data.T(i1,:) = prod(squeeze(img.get('T')))^0.5;
%     r.data.pos(i1,:) = img.get('pos');

    figure(98);clf;
    subplot(1,2,1);
    plot(r.data.param(1:i1),r.data.N,'o');
    ylim([0,Inf]);
    grid on
    subplot(1,2,2);
%     plot(r.data.param(1:i1),r.data.T,'o');
%     ylim([0,250e-9]);
    plot(r.data.param(1:i1),r.data.N(:,1)./sum(r.data.N,2),'o-');
    grid on
%     subplot(1,2,3);
%     plot(r.data.param(1:i1),r.data.becFrac,'o');
%     ylim([0,Inf]);
%     grid on
%     subplot(1,2,4);
%     plot(r.data.param(1:i1),r.data.becFrac,'o');
%     ylim([0,Inf]);
%     grid on
%     subplot(1,3,3);
%     peak_density = r.data.N.*((0.5*const.muB*110e-2)./(2*const.kb*r.data.T)).^3;
%     average_speed = sqrt(16*const.kb*r.data.T/(pi*const.mRb));
%     cross_section = 8*pi*(100*const.aBohr)^2;
%     dbwavelength = sqrt(2*pi*const.hbar^2./(const.mRb*const.kb*r.data.T));
%     r.data.psd = peak_density.*dbwavelength.^3;
%     r.data.col_rate = peak_density.*average_speed.*cross_section/8;
%     plot(r.data.param(1:i1),r.data.col_rate,'o');
%     grid on;
end


end