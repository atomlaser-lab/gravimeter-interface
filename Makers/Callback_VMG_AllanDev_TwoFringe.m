function Callback_VMG_AllanDev_TwoFringe(r)
if r.isInit()
    r.data.T =  [1:1:10000]; %for dipoles: [1.5:0.5:8] %for RF: [0.5:0.25:2]
    r.data.phase = [10e-6,40e-6];% phase AOM %19
    r.c.setup('var',r.data.phase,r.data.T);
elseif r.isSet()
    r.make(r.devices.opt,'param1',r.data.phase(r.c(1)),'param2',r.data.T(r.c(2))).upload;
    fprintf(1,'Run %d/%d, Phase = %.3f Degrees\n',r.c.now,r.c.total,2*r.data.phase(r.c(1)));


elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis_DualState('last',1);

    %plotting
    r.data.files{i1,1} = img(1).raw.files;
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');
    r.data.peakOD(i1,:) = img.get('peakOD');
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:),2);
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:),2);

    fid = fopen('D:\data\VMG_autosave\raman beam positions\poynting_vector_binary_indictor.txt','w'); fprintf(fid,'1'); fclose(fid);    
    pause(1);
    img_no = sscanf(r.data.files{i1}.name,'bec%d.bin');
    imgHyGG = sprintf('D:\\data\\VMG_autosave\\raman beam positions\\%g_HyGG.tif',img_no);
    imgG = sprintf('D:\\data\\VMG_autosave\\raman beam positions\\%g_G.tif',img_no);
    img_DRK = sprintf('D:\\data\\VMG_autosave\\raman beam positions\\%g_DRK.tif',img_no);
    temp = double(imread(imgHyGG) - imread(img_DRK));
    [tempx,tempy] = find(max(max(temp)) == temp);
    r.data.HyGGx_spot(i1) = tempx(1)*5.5e-6;
    r.data.HyGGy_spot(i1) = tempy(1)*5.5e-6;
    temp = double(imread(imgG) - imread(img_DRK));
    [f,g] = Raman_Gaussian_Fit(1:512,1:500,temp);
    r.data.Gx_spot(i1) = f.b*5.5e-6;
    r.data.Gy_spot(i1) = f.d*5.5e-6;

%     figure(133);clf;
%     plot(1:i1,r.data.Rsum(:,1)-r.data.Rsum(:,2),'.');
%     enhformat('parameter [a.u]','Population Difference')
%     grid on;
%     ylim([-1 1]);

    if r.c.done(1)
        %saving to VMG_autosave folder in D drive
        data = r.data;
        save('D:\data\VMG_autosave\data.mat','data');
        saveas(gcf,'D:\data\VMG_autosave\data.fig');
    end
end