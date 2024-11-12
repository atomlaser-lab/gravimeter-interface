function Callback_VMG_AllanDev_TwoFringe(r)
if r.isInit()
    r.data.T =  [1:1:10000]; %for dipoles: [1.5:0.5:8] %for RF: [0.5:0.25:2]
    r.data.phase = [13.2979,103.2979];% phase AOM %19
    r.c.setup('var',r.data.phase,r.data.T);
elseif r.isSet()
    r.make(r.devices.opt,'param1',r.data.phase(r.c(1)),'param2',r.data.T(r.c(2))).upload;
    fprintf(1,'Run %d/%d, Phase = %.3f Degrees\n',r.c.now,r.c.total,2*r.data.phase(r.c(1)));
    if r.c.done(1)
        %saving to VMG_autosave folder in D drive
        data = r.data;
        save('D:\data\VMG_autosave\data.mat','data');
        saveas(gcf,'D:\data\VMG_autosave\data.fig');
    end
end