function sq = trigger_ndi(r,opt)


sq = initSequence;
sq.delay(1);
sq.camDelay = sq.time;
sq.delay(1);
makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_time,'cycle time',500e-3,...
    'imaging amplitude',opt.nd.pulse_power,'num_images',opt.nd.ref_images,...
    'pulse delay',opt.nd.pulse_delay);

sq.delay(1);
makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_time,'cycle time',opt.nd.cycle_time,...
    'imaging amplitude',opt.nd.pulse_power,'num_images',opt.nd.num_images(1),...
    'pulse delay',opt.nd.pulse_delay);

sq.delay(1);
makeImagingSequence(sq,'tof',opt.tof,'pulse time',40e-6,'repump delay',100e-6,...
    'repump time',200e-6,'cam time',5e-6,'cycle time',100e-3,...
    'manifold',1,'imaging freq',0,'imaging amplitude',0.5,...
    'repump shutter delay',2e-3,'imaging_field',1,'image type','horizontal');

setSafeValues(sq);

r.sq = sq;
r.urun;

end