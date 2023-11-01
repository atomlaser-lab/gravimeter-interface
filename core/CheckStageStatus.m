function CheckStageStatus(DesiredEndStage)

if strcmpi(DesiredEndStage,'CMOT') == 1
    opt.CMOT_status = 1;
elseif strcmpi(DesiredEndStage,'PGC') == 1
    opt.CMOT_status = 1;
    opt.PGC_status = 1;
elseif strcmpi(DesiredEndStage,'MagLoad') == 1
    opt.CMOT_status = 1;
    opt.PGC_status = 1;
    opt.LoadMagTrap_status = 1;
elseif strcmpi(DesiredEndStage,'MagEvap') == 1
    opt.CMOT_status = 1;
    opt.PGC_status = 1;
    opt.LoadMagTrap_status = 1;
    opt.MagEvaporation_status = 1;
elseif strcmpi(DesiredEndStage,'OpticalLoad') == 1
    opt.CMOT_status = 1;
    opt.PGC_status = 1;
    opt.LoadMagTrap_status = 1;
    opt.MagEvaporation_status = 1;
    opt.LoadOpticalTrap_status = 1;
elseif strcmpi(DesiredEndStage,'OpticalEvap') == 1
    opt.CMOT_status = 1;
    opt.PGC_status = 1;
    opt.LoadMagTrap_status = 1;
    opt.MagEvaporation_status = 1;
    opt.LoadOpticalTrap_status = 1;
    opt.OpticalEvaporation_status = 1;
end


if opt.CMOT_status == 0
    opt.PGC_status = 0;
    opt.LoadMagTrap_status = 0;
    opt.MagEvaporation_status = 0;
    opt.LoadOpticalTrap_status = 0;
    opt.OpticalEvaporation_status = 0;
    opt.BECCompression_status = 0;
    opt.MagneticInsensitive_status = 0;
elseif opt.PGC_status == 0
    opt.LoadMagTrap_status = 0;
    opt.MagEvaporation_status = 0;
    opt.LoadOpticalTrap_status = 0;
    opt.OpticalEvaporation_status = 0;
    opt.BECCompression_status = 0;
    opt.MagneticInsensitive_status = 0;
elseif opt.LoadMagTrap_status == 0
    opt.MagEvaporation_status = 0;
    opt.LoadOpticalTrap_status = 0;
    opt.OpticalEvaporation_status = 0;
    opt.BECCompression_status = 0;
    opt.MagneticInsensitive_status = 0;
elseif opt.MagEvaporation_status == 0
    opt.LoadOpticalTrap_status = 0;
    opt.OpticalEvaporation_status = 0;
    opt.BECCompression_status = 0;
    opt.MagneticInsensitive_status = 0;
elseif opt.LoadOpticalTrap_status == 0
    opt.OpticalEvaporation_status = 0;
    opt.BECCompression_status = 0;
    opt.MagneticInsensitive_status = 0;
elseif opt.OpticalEvaporation_status == 0
    opt.OpticalEvaporation_status = 0;
    opt.BECCompression_status = 0;
    opt.MagneticInsensitive_status = 0;
end
end