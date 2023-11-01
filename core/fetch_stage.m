function cooling_stage = fetch_stage()
    % Fetch 'opt' structure from the base workspace
    opt = evalin('base', 'opt');

    % List of fields to check in the order given
    fields = {
        'MOT_status', 'MOT';
        'CMOT_status', 'CMOT';
        'PGC_status', 'PGC';
        'LoadMagTrap_status', 'Load in the mag Trap';
        'MagEvaporation_status', 'Mag evaporation';
        'LoadOpticalTrap_status', 'Load in the optical trap';
        'OpticalEvaporation_status', 'Optical evaporation'
    };

    % Start with the assumption that no stages are complete
    cooling_stage = 'None';

    % Iterate over the fields and check their status
    for i = 1:size(fields, 1)
        if opt.(fields{i, 1}) == 1
            cooling_stage = fields{i, 2};
        else
            % Stop checking when a field is found to be 0
            break;
        end
    end
end
