function V = TrapPtoV(light_type,relative_power)

if relative_power > 1
    relative_power = 1;
elseif relative_power < 0
    relative_power = 0;
end

Vi = linspace(0,10,1e2);

if strcmpi(light_type,'trap')
    data = [8	105;
            7	99;
            6	81;
            5	54;
            5.50000000000000	67;
            4.50000000000000	41;
            3.50000000000000	17;
            4	27.5000000000000;
            3	8.30000000000000;
            2.50000000000000	2.70000000000000;
            2	0.360000000000000;
            2.10000000000000	0.610000000000000;
            2.20000000000000	0.950000000000000;
            2.30000000000000	1.40000000000000;
            2.40000000000000	1.94000000000000;
            2.60000000000000	3.40000000000000;
            2.70000000000000	4.40000000000000;
            2.80000000000000	5.50000000000000;
            2.90000000000000	6.60000000000000;
            1.50000000000000	0];
    data(:,2) = data(:,2)./max(data(:,2));
    V = interp1(data(:,2),data(:,1),relative_power,'pchip');
    
elseif strcmpi(light_type,'repump')
    data = [6 3.9;
            5 3.8;
            4 3.3;
            3 2.0;
            2.5 1.15;
            2 0.43;
            1.7 0.13;
            1.5 0.03];
    data(:,2) = data(:,2)./max(data(:,2));
    V = interp1(data(:,2),data(:,1),relative_power,'pchip');
elseif strcmpi(light_type,'nd')
    data = [8.0000    1.5000;
            7.0000    1.4570;
            6.0000    1.1930;
            5.0000    0.8100;
            4.0000    0.4200;
            3.0000    0.1190;
            2.0000    0.0046;
            3.5000    0.2450;
            4.5000    0.6090;
            5.5000    1.0000];
    [~,k] = sort(data(:,1));
    data = data(k,:);
    data(:,2) = data(:,2)/max(data(:,2));
    V = interp1(data(:,2),data(:,1),relative_power,'pchip');
    V(relative_power == 0) = 0;
elseif strcmpi(light_type,'image')
    data = [8,1.24;
            7,1.39;
            6,1.42;
            5,1.18;
            4.5,0.96;
            4,0.72;
            3.5,0.47;
            3,0.25;
            2.75,0.16;
            2.5,0.1;
            2.25,0.052;
            2,0.026;
            1.5,0.014;
            0,0.014];
    data = data - data(end,2);
    data = data(data(:,1) <= 6 & data(:,1) > 0,:);
    data(:,2) = data(:,2)./max(data(:,2));
    V = interp1(data(:,2),data(:,1),relative_power,'pchip');
    V(relative_power == 0) = 0;
end

end
