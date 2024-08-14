function Callback_StabilityAnalysis(r)

if r.isInit()
%     r.data.power = 0.328;
    
    r.c.setup(Inf);
elseif r.isSet()
    
%     r.make(25.5,730e-3,1.48,r.data.power,0,130e-3);
%     r.upload;
    fprintf(1,'Run %d/%d\n',r.c.now,r.c.total);
%     pause(10);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1);
%     img = Abs_Analysis('last');
%     if ~img(1).raw.status.ok()
%         %
%         % Checks for an error in loading the files (caused by a missed
%         % image) and reruns the last sequence
%         %
%         r.c.decrement;
%         return;
%     end
    %
    % Store raw data
    %
%     r.data.img{i1,1} = img;
%     r.data.files{i1,1} = img.raw.files;
%     r.data.date(i1,1) = img.raw.files(1).datenum;
    %
    % Get processed data
    %
%     r.data.N(i1,:) = img.get('N');
%     r.data.Nsum(i1,:) = img.get('Nsum');
%     r.data.width(i1,:) = img.clouds(1).becWidth;
%     r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
%     r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));

    [~,N,dout] = FMI_Analysis;
    r.data.N(i1,:) = [N.N1,N.N2];
    r.data.R(i1,1) = N.R;
    r.data.Nsum(i1,:) = N.sum;
    r.data.Rsum(i1,1) = N.Rsum;
    r.data.date(i1,1) = now;
    r.data.dout(i1,1) = dout;
    
%     r.data.fmi.t(:,i1) = nlf.x;
%     r.data.fmi.v(:,i1) = nlf.y;
%     
%     if isfield(r.devices,'d2')
%         d = r.devices.d2.getRAM;
%         nlf = nonlinfit(d.t,d.data(:,1));
%         nlf.useErr = false;
%         gauss = @(A,w,x0,x) A*exp(-(x-x0).^2/w^2);
%         nlf.setFitFunc(@(At,C,w,x0,f0,phi,y0,x) gauss(At,w,x0,x).*(1 + C*sin(2*pi*f0*x+phi)) + y0);
%         nlf.bounds2('At',[0,1,max(nlf.y)/2],'C',[0,1,0.5],'w',[5e-6,150e-6,10e-6],...
%             'x0',[2e-4,3e-4,2.5e-4],'f0',[0.95e6,1.1e6,1.04e6],'phi',[-5*pi,5*pi,0],...
%             'y0',[-0.05,0.05,nlf.y(1)]);
%         nlf.fit;
%         figure(1);clf;
%         plot(nlf.x,nlf.y,'.-');
%         hold on
%         plot(nlf.x,nlf.f(nlf.x),'-');
% 
%         r.data.c(:,i1) = nlf.c(:,1);
%         r.data.pd.t(:,i1) = nlf.x;
%         r.data.pd.v(:,i1) = nlf.y;
%     end

%     d = DataAcquisition('192.168.1.21');
%     d.fetch;
%     d.getRAM;
%     r.data.t{i1,1} = d.t;
%     r.data.v{i1,1} = d.data;
    
    figure(98);clf;
%     subplot(1,2,1);
    plot((r.data.date(1:i1) - r.data.date(1))*24*60,r.data.R(1:i1),'o-');
    hold on
    plot((r.data.date(1:i1) - r.data.date(1))*24*60,r.data.Rsum(1:i1),'sq-');
    hold off;
    plot_format('Time since start [min]','Number','',12);
    grid on;
%     ylim([0,1]);


    
    
end