function Callback_MeasureImagingFrequency(r)

if r.isInit()
    r.data.detuning = const.randomize([-16:2:16,-5:2:5]);
%     r.data.detuning = const.randomize(repmat([-20:1:20,-6.5:1:6.5],4,1));
    r.data.param = 5;
    r.c.setup('var',r.data.detuning,r.data.param);
elseif r.isSet()
    r.make(r.devices.opt,'detuning',r.data.detuning(r.c(1))).upload;
    fprintf(1,'Run %d/%d, Detuning = %.3f MHz, Param = %.3f\n',r.c.now,r.c.total,r.data.detuning(r.c(1)),r.data.param(r.c(2)));
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis('last',1);
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
    
    r.data.files{i1,i2} = img.raw.files;
    r.data.N(i1,i2) = img.get('N');
    r.data.OD(i1,i2) = img.get('peakOD');

    r.data.files{i1,1} = img.raw.files;
    r.data.N(i1,1) = img.get('N');
    r.data.pos(i1,:) = squeeze(img.get('pos'));
    figure(123);clf;
    plot(r.data.detuning(1:i1),r.data.N,'o');
    plot_format('Imaging detuning [MHz]','Number [arb units]','',12);
    grid on;
    ylim([0,Inf]);

    if r.c.done(1) || i1 > 6
%         nlf = nonlinfit(r.data.detuning(1:i1),r.data.N(1:i1)/1e6,0.05);
%         nlf.setFitFunc(@(A1,w1,x1,A2,w2,x2,x) A1./(1 + 4*(x-x1).^2/w1^2) + A2./(1 + 4*(x-x2).^2/w2^2));
%         [~,idx] = max(nlf.y);
%         nlf.bounds2('A1',[0,1e3,max(nlf.y)],'w1',[1,12,6],'x1',[min(nlf.x),max(nlf.x),nlf.x(idx)],...
%             'A2',[0,1e3,max(nlf.y)],'w2',[1,12,6],'x2',[min(nlf.x),max(nlf.x),nlf.x(idx) - 5]);
% % %         nlf.setFitFunc(@(A1,w1,x1,A2,x2,x) A1./(1 + 4*(x-x1).^2/w1^2) + A2./(1 + 4*(x-x2).^2/w1^2));
% % %         [~,idx] = max(nlf.y);
% % %         nlf.bounds2('A1',[0,1e3,max(nlf.y)],'w1',[1,12,6],'x1',[min(nlf.x),max(nlf.x),nlf.x(idx)],...
% % %             'A2',[0,1e3,max(nlf.y)],'x2',[min(nlf.x),max(nlf.x),nlf.x(idx) - 5]);
%         nlf.ex = (nlf.y == 0);
%         nlf.fit;
%         hold on;
%         xplot = linspace(min(nlf.x),max(nlf.x),1e2);
%         plot(xplot,nlf.f(xplot)*1e6,'--','linewidth',2);
%         r.data.nlf = nlf;
%         s1 = sprintf('A1 = %.1f, x1 = %.1f MHz, w1 = %.1f MHz\n',nlf.get('A1',1),nlf.get('x1',1),nlf.get('w1',1));
%         s2 = sprintf('A2 = %.1f, x2 = %.1f MHz, w2 = %.1f MHz\n',nlf.get('A2',1),nlf.get('x2',1),nlf.get('w2',1));
%         xmax = fminbnd(@(x) -nlf.f(x),min(nlf.x),max(nlf.x));
%         fmax = nlf.f(xmax);
%         s3 = sprintf('Max response is %.3f at %.3f MHz',fmax,xmax);
%         text(0.05,0.9,[s1,s2,s3],'units','normalized');

        nlf = nonlinfit(r.data.detuning(1:i1),r.data.N(1:i1)/1e6,0.05);
        nlf.ex = (nlf.y == 0);% | (nlf.x == 10);
        nlf.setFitFunc(@(A,w,x0,x) A./(1 + 4*(x - x0).^2/w^2));
        nlf.bounds2('A',[0,2*max(nlf.y),max(nlf.y)],'w',[0,10,6],'x0',[-3,3,0]);
        nlf.fit;
        hold on
        xplot = linspace(min(nlf.x),max(nlf.x),1e2);
        plot(xplot,nlf.f(xplot)*1e6,'--','linewidth',2);
        title(sprintf('Resonance = %.1f MHz, Width = %.1f MHz',nlf.get('x0',1),nlf.get('w',1)));
        r.data.nlf = nlf;
    end
end


end