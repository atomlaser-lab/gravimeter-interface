function NeededVoltage = DipolePtoV(DipoleType,DesiredPower)

% Powers measured on 26/04/2022
% RedPower has P_max = 20 W
% Raycus has P_max = 16 W at 3.1 V


if ~(strcmpi(DipoleType,'RedPower') || strcmpi(DipoleType,'Raycus') || strcmpi(DipoleType,'DKC') || strcmpi(DipoleType,'feedback'))
    error('Must be RedPower or Raycus')
end

if strcmpi(DipoleType,'RedPower')
    if any(DesiredPower > 20 | DesiredPower < 0)
        error('Power cannot be out of range [0,20] W!');
    end
    NeededVoltage = (DesiredPower + 1.407)/2.5341; %as measured on 19/11/2024

    if NeededVoltage <= 0
        NeededVoltage = 0;
    elseif NeededVoltage > 10
        NeededVoltage = 10;
    end

elseif strcmpi(DipoleType,'Raycus')
    %
    % Linear regression done on data taken on 19/11/2024
    %
    if DesiredPower > 14.1
        error('Power must be less than 14.1 W');
    elseif DesiredPower < 0
        error('Power must be larger than 0 W');
    else
        NeededVoltage = (DesiredPower + 1.3069)/4.4077; %as measured on 19/11/2024
    end

elseif strcmpi(DipoleType,'DKC')
    if DesiredPower > 11.5
        error('Power must be less than 11.5 W');
    elseif DesiredPower < 0
        error('Power must be larger than 0 W');
    else
        x = DesiredPower;
        NeededVoltage = (x + 0.9196)/3.18;
        if NeededVoltage <= 0
            NeededVoltage = 0;
        end
    end

elseif strcmpi(DipoleType,'Feedback')
    if DesiredPower > 12
        error('Power must be less than 1.2 W');
    elseif DesiredPower < 0
        error('Power must be larger than 0 W');
    else
        x = DesiredPower;
        NeededVoltage = (x + 1.23)/3.73;
        if NeededVoltage < 0
            NeededVoltage = 0;
        end
    end
end


end




