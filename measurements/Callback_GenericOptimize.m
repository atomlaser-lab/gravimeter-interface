function Callback_GenericOptimize(r)
xAxis = 'Load (s)';
title = 'MOT Load Rate';

if r.isInit()
    %Initialize run

    %scan these variables (scans through param then param2)
    r.data.freq = 10:5:35;
    r.data.power = [2:0.2:4,4.5,5,6];

    r.c.setup('var',r.data.power,r.data.freq);

%     if r.data.counter == 1
%         r.c.setup('var',r.data.param);
%     elseif r.data.counter == 2
%         r.c.setup('var',r.data.param,r.data.param2); %note the order here determines what varargin{i} is
%     end

elseif r.isSet()
    %print parameters after each run
    if r.data.counter == 1
        r.make(r.data.param(r.c(1)));
        r.upload;
        fprintf(2,'Run %d/%d, Param = %.3f\n',r.c.now,r.c.total,...
            r.data.param(r.c(1)));

    elseif r.data.counter == 2
        r.make(r.data.param(r.c(1)),r.data.param2(r.c(2)));
        r.upload;
        fprintf(2,'Run %d/%d, Param = %.3f, Param2 = %.3f\n',r.c.now,r.c.total,...
            r.data.param(r.c(1)),r.data.param2(r.c(2)));
    end
elseif r.isAnalyze()

    if r.data.counter == 1
        i1 = r.c(1);
        img = Abs_Analysis('last');


        r.data.files{i1} = img.raw.files;
        r.data.N(i1) = img.get('N');
        r.data.T(i1) = sqrt(prod(img.clouds.T));
        r.data.OD(i1) = img.clouds.peakOD;
        r.data.GaussWidth(i1,:) = img.get('gaussWidth');
        r.data.pos(i1,:) = squeeze(img.get('pos'));

    elseif r.data.counter == 2

        img = Abs_Analysis('last');
        i1 = r.c(1);
        i2 = r.c(2);

        r.data.files{i1,i2,:} = img.raw.files;
        r.data.N(i1,i2,:) = img.get('N');
        r.data.T(i1,i2,:) = sqrt(prod(img.clouds.T));
        r.data.OD(i1,i2,:) = img.clouds.peakOD;
        r.data.GaussWidth(i1,i2,:) = img.get('gaussWidth');
    end

    %make plot during run for Varargin = 1
    if r.data.counter == 1
        %position
%         figure(2)
%         plot(r.data.param(1:1:size(r.data.N,2)),r.data.pos,'o')
%         xlabel(xAxis)
%         ylabel('Pos')
        %Number
        figure(3)
        subplot(2,1,1)
        plot(r.data.param(1:1:size(r.data.N,2)),r.data.N)
        xlabel(xAxis)
        ylabel('Fitted Atom Number')
        subplot(2,1,2)
        plot(r.data.param(1:1:size(r.data.N,2)),r.data.OD)
        xlabel(xAxis)
        ylabel('Max OD')
        sgtitle(title)
        %Width
%         figure(4)
%         plot(r.data.param(1:1:size(r.data.N,2)),r.data.GaussWidth)
%         xlabel(xAxis)
%         ylabel('Fitted Cloud Size')
%         sgtitle(title)

    end




    if r.data.counter == 2

        %Get Temperature for each param2
        if r.c(1) == size(r.data.param,2)
            for nn = 1:2
                w = squeeze(r.data.GaussWidth(:,r.c(2),nn));
                lf = linfit(r.data.param,w.^2,2*w.*20e-6);
                lf.setFitFunc('poly',[0,2]);
                lf.fit;
                figure(5);clf
                lf.plot;
                r.data.Tfit(i2,nn) = lf.c(2,1)*const.mRb/const.kb*1e6;
                r.data.Terr(i2,nn) = lf.c(2,2)*const.mRb/const.kb*1e6;
            end
        end

        %Plot temperature vs param2
        if r.c(2) >= 2
            figure(6);clf
            plot(r.data.param2(1:1:size(r.data.Tfit(:,1))),r.data.Tfit(:,1))
            hold on
            plot(r.data.param2(1:1:size(r.data.Tfit(:,1))),r.data.Tfit(:,2))
            ylabel('Temperature (uK)')
            xlabel(xAxis)
            legend('X Temperature','Y Temperature')
            sgtitle(title)
        end
        
        if r.c(1) == size(r.data.param,2)
%             r.data.N(i1,i2,:) = img.get('N');
%             plot(r.data.param2(1:1:size(r.data.N(4,:),2)),r.data.N(4,:))
        end

    end


end