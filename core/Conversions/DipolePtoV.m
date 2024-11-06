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

    %%% OLD data
    %     data = [0.882	0.910;
    %             0.803	0.710;
    %             0.603	0.230;
    %             0.673	0.390;
    %             0.713	0.480;
    %             1.013	1.250;
    %             0.990	1.190;
    %             0.950	1.090;
    %             0.900	0.960;
    %             1.300	2.000;
    %             1.190	1.730;
    %             1.150	1.620;
    %             1.440	2.370;
    %             1.670	3.000;
    %             1.600	2.810;
    %             1.520	2.600;
    %             2.000	3.870;
    %             2.500	5.150;
    %             2.100	4.130;
    %             2.200	4.390;
    %             1.900	3.610;
    %             1.950	3.730;
    %             2.600	5.460;
    %             2.750	5.880;
    %             3.000	6.590;
    %             3.500	7.960;
    %             4.000	9.270;
    %             4.500	10.60;
    %             5.000	11.90;
    %             5.500	13.20;
    %             6.000	14.60;
    %             6.500	15.90;
    %             7.000	17.20;
    %             7.500	18.40;
    %             8.000	19.60;
    %             8.500	20.70];
    %
    %    data = [0.8	0.7
    %         0.9	0.99
    %         1	1.26
    %         1.1	1.54
    %         1.2	1.82
    %         1.3	2.09
    %         1.4	2.37
    %         1.5	2.66
    %         1.6	2.95
    %         1.7	3.23
    %         1.8	3.5
    %         1.9	3.78
    %         2	4.05
    %         2.3	4.87
    %         2.6	5.72
    %         2.9	6.61
    %         3.2	7.47
    %         3.5	8.33
    %         3.8	9.17
    %         4.1	9.99
    %         4.4	10.8
    %         4.7	11.6
    %         5	12.4
    %         5.3	13.2
    %         5.6	14.1
    %         5.9	14.9
    %         6.2	15.7
    %         6.5	16.5
    %         6.8	17.4
    %         7.1	18.1
    %         7.3	18.7
    %         0.7	0.43
    %         0.6	0.18
    %         1.05	1.38
    %         1.15	1.66
    %         1.25	1.94
    %         1.35	2.22
    %         1.45	2.5]; % COLLECTED 2024-04-11

    % NeededVoltage = interp1(data(2:end,2),data(2:end,1),DesiredPower,'pchip');

    %
    % Linear regression done on data taken on 03/11/2023
    %
    % NeededVoltage = (DesiredPower + 1.4697)/2.7425;
    %     NeededVoltage = (DesiredPower + 1.4990)/2.778;
    NeededVoltage = (DesiredPower + 1.5494)/2.7587; %as measured on 23/05/24

    %
    % Cubic fit, closed-loop control 24/04/2024
    %
    % NeededVoltage = -0.00822 + 0.1787*DesiredPower - 0.00143782*DesiredPower.^2 - 1.2626e-5*DesiredPower.^3;

    %     %%% Linear fit:
    %     NeededVoltage = (DesiredPower + 0.286)/6.348;


    if NeededVoltage < 0
        NeededVoltage = -0.1;
    elseif NeededVoltage > 10
        NeededVoltage = 10;
    elseif DesiredPower > 17
        error('Closed-loop RedPower control is limited to 17 W!');
    elseif DesiredPower == 0
        NeededVoltage = -0.1;
    end

elseif strcmpi(DipoleType,'Raycus')
    % %%% OLD data
    %     data = [ 0.00000,  0.40000;
    %              0.25000,  0.40000;
    %              0.50000,  1.18000;
    %              0.75000,  2.24000;
    %              1.00000,  3.23000;
    %              1.50000,  5.00000;
    %              2.00000,  5.84000;
    %              2.50000,  6.39000;
    %              2.90000,  8.73000;
    %              2.25000,  5.97000;
    %              1.25000,  4.16000;
    %              1.75000,  5.62000;
    %              2.75000,  7.25000;
    %              3.00000,  9.72000;
    %              3.20000,  11.00000;
    %              3.50000,  12.60000];
    %
    %     % New data for Raycus, as of 28/08/2023
    %     data = [0.3	0.115;
    %             0.35 0.336;
    %             0.4	 0.559;
    %             0.5	 0.995;
    %             1	 3.257;
    %             2	 7.747;
    %             3	 12.097;
    %             3.9	 16.297];
    %
    % New data for Raycus, as of 01/09/2023
    %     data = [0.3	0.2065/2;
    %             0.35 0.5105/2;
    %             0.4	 0.8165/2;
    %             0.5	 1.4745/2;
    %             1	 4.5045/2;
    %             2	 10.3945/2;
    %             3	 15.5945/2;
    %             3.1	 16.0945/2];
    %
    %     if DesiredPower < 0.4
    %         NeededVoltage = 0;
    %     elseif DesiredPower > max(data(:,2))
    %         error('Power must be less than %.3f!',max(data(:,2)));
    %     else
    %         NeededVoltage = interp1(data(2:end,2),data(2:end,1),DesiredPower,'pchip');
    %     end

    %
    % Linear regression done on data taken on 03/11/2023
    %
    if DesiredPower > 11
        error('Power must be less than 11 W');
    elseif DesiredPower < 0
        error('Power must be larger than 0 W');
    else
        % NeededVoltage = (DesiredPower + 1.81)./6.7235;    %%% OLD
        %         NeededVoltage = (DesiredPower + 1.4350)./6.1529;
        %         NeededVoltage = (DesiredPower - 0.15)/3.33; %Closed-loop 26/04/2024
        NeededVoltage = (DesiredPower + 1.458)/5.8594; %as measured on 23/05/24
    end

elseif strcmpi(DipoleType,'DKC')
    if DesiredPower > 11.5
        error('Power must be less than 11.5 W');
    elseif DesiredPower < 0
        error('Power must be larger than 0 W');
    else
        x = DesiredPower;
        NeededVoltage = (x + 0.9196)/3.18;
        if NeededVoltage < 0
            NeededVoltage = 0;
        end
    end

elseif strcmpi(DipoleType,'Feedback')
    if DesiredPower > 1.5
        error('Power must be less than 1.5 W');
    elseif DesiredPower < 0
        error('Power must be larger than 0 W');
    else
        x = DesiredPower;
        % NeededVoltage = (x + 0.338)/1.1562;
        NeededVoltage = (x + 0.6273)/2.1042;
        if NeededVoltage < 0
            NeededVoltage = 0;
        end
    end
end


end




