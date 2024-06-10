%% yosri configuration
function Yosri ()
% Modify existing variables in the workspace
    global opt m r const

opt.tof = 217.5e-3;
opt.MOT_LoadTime = 12;
opt.TwoStateImaging = 0;
% f1 = (const.f_Rb_groundHFS/1e6 - 315e-3 + (17-5.3)*1e-3)/2*1e6;
% m.writeList(f1,(const.f_Rb_groundHFS)/2);
r.makerCallback = @Alternative_NewSF_300524;
% Update the modified objects in the workspace
    assignin('base', 'opt', opt);
    assignin('base', 'm', m);
    assignin('base', 'r', r);
    
end