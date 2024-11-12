function Callback_Vis_vs_T(r)

if r.isInit()
    r.data.T =  [500e-6 750e-6 1000e-6]; %for dipoles: [1.5:0.5:8] %for RF: [0.5:0.25:2]
    r.data.phase = const.randomize(linspace(0,180,181));% phase AOM %19
    r.data.limit = 1000;
    r.c.setup('var',r.data.phase,r.data.T);
elseif r.isSet()
    r.make(r.devices.opt,'param1',r.data.phase(r.c(1)),'param2',r.data.T(r.c(2))).upload;
    fprintf(1,'Run %d/%d, Phase = %.3f Degrees, T = %.1f ms\n',r.c.now,r.c.total,2*r.data.phase(r.c(1)),1e3*r.data.T(r.c(2)));
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis_DualState('last',1);

    r.data.files{i1,i2,1} = img(1).raw.files;
    r.data.N(i1,i2,:) = img.get('N');
    r.data.Nsum(i1,i2,:) = img.get('Nsum');
    r.data.peakOD(i1,i2,:) = img.get('peakOD');
    r.data.R(i1,i2,:) = r.data.N(i1,i2,:)./sum(r.data.N(i1,i2,:));
    r.data.Rsum(i1,i2,:) = r.data.Nsum(i1,i2,:)./sum(r.data.Nsum(i1,i2,:));

if (~img(1).raw.status.ok() || r.data.N(i1,1) == 0 || r.data.N(i1,1) > 1e7 || r.data.Rsum(i1,1) > 0.99 || (r.data.N(i1,1) == 0 && r.data.N(i1,2) == 0))
        %         if (~img(1).raw.status.ok() || r.data.N(i1,1) == 0)
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        if (~img(1).raw.status.ok())
            warning('Imaging failed!')
            r.c.decrement;

        elseif ( r.data.N(i1,1) == 0 || r.data.N(i1,1) > 1e7 || r.data.Rsum(i1,1) > 0.99 || (r.data.N(i1,1) == 0 && r.data.N(i1,2) == 0))
            %             elseif (r.data.N(i1,1) == 0)
            r.data.limit=1+r.data.limit;
            warning('Number of atoms too low!')
            r.c.decrement;
            if r.data.limit == 100
                r.data.limit = 0;
                r.stop
                error('Run failed! lock your lasers.')
            end
        end
        return;
    elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end

    figure(42)
    clf
    hold on
    scatter(2*r.data.phase(1:i1),r.data.Rsum([1:i1],i2,1) - r.data.Rsum([1:i1],i2,2),'filled')
    xlim([0 360])
    ylim([-1 1])
    xlabel('LG Rotation Phase (Deg)')
    ylabel('N_1-N_2')
    title(sprintf('T = %.2f',1e6*r.data.T(i2)))

    if r.c(1) == length(r.data.phase)
%         [f,g] = createFit(r.data.x',(0.5*(1+y1))');
%         x=linspace(0,360);
%         y = f.a*sin(pi*x/180+f.b)+f.c;
%         plot(x,y,'Linewidth',2,'Color',"blue");
        saveas(gcf,sprintf('D:\\data\\VMG_autosave\\T =  %.3f ms, OAM = 5000.png',1e3*r.data.T(i2)));
%         r.data.vis(i2) = 2*f.a;
% 
%         figure(43);
%         plot(r.data.raycus(1:i2),r.data.vis,'o')
%         ylim([0 1])
%         grid on;
%         enhformat('','Visibility')
    end

    if r.c.done(1)

        %saving to VMG_autosave folder in D drive
        data = r.data;
        save('D:\data\VMG_autosave\data.mat','data');
        saveas(gcf,'D:\data\VMG_autosave\data.fig');

    end

end