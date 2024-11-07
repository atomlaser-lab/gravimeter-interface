function [out] = FtoV(LightType,Frequency)
% Interpolation range
Voltage = (0:0.01:10);

% Voltage to Frequency Functions
% TrappingFrequency = 2*(53.051+8.6164*(Voltage)-1.5183*((Voltage).^2)+.24203*((Voltage).^3)-.010976*Voltage.^4)-211.79;
% % RepumpFrequency = 2*(51.919+8.5694*(Voltage)-1.5263*(Voltage.^2)+.24217*(Voltage.^3)-.010947*Voltage.^4) -211.79;
% RepumpFrequency = 55.1531 + 5.514*Voltage - 0.3122*Voltage.^2 + 0.048*Voltage.^3 - 156.947/2;
% % ImagingFrequency = 0.5168+16.0185*(Voltage-8.386)-0.112*(Voltage-8.386).^2;
% ImagingFrequency = 2*(50.3705 + Voltage*5.2047 + 0.1675*Voltage.^2) - 0.5*(266.65 + 156.947);

%
% Measurements from 2024-11-08
%
RepumpFrequency = (51.3531 + 9.5677*Voltage - 1.8097*Voltage.^2 + 0.2743*Voltage.^3 - 0.0121*Voltage.^4) - 0.5*156.946;
ImagingFrequency = 2*(51.3753 + 9.3530*Voltage - 1.7427*Voltage.^2 + 0.2667*Voltage.^3 - 0.0119*Voltage.^4) - 0.5*(266.65 + 156.947);
TrappingFrequency = 2*(53.5097 + 8.1718*Voltage - 1.3813*Voltage.^2 + 0.2251*Voltage.^3 - 0.0102*Voltage.^4) - 0.5*(266.65 + 156.947);
Trapping2DFrequency = 2*(52.5931 + 8.5238*Voltage - 1.5314*Voltage.^2 + 0.2462*Voltage.^3 - 0.0112*Voltage.^4) - 0.5*(266.65 + 156.947);
PushFrequency = 2*(51.3753 + 9.3530*Voltage - 1.7427*Voltage.^2 + 0.2667*Voltage.^3 - 0.0119*Voltage.^4) - 0.5*(266.65 + 156.947);
%
% use input string, the input desired frequency and the above functions to get the required voltage
%
switch lower(LightType)
    case 'trap'
        out = interp1(TrappingFrequency,Voltage,-Frequency);
    case 'repump'
        out = interp1(RepumpFrequency,Voltage,-Frequency);
    case {'image','imaging'}
        out = interp1(ImagingFrequency,Voltage,-Frequency);
    case '2d'
        out = interp1(Trapping2DFrequency,Voltage,-Frequency);
    case 'push'
        out = interp1(PushFrequency,Voltage,-Frequency);
    otherwise
        error('First argument must be one of ''trap'', ''repump'', ''image'', ''2d'', or ''push''');     
end

% 
% Check if frequency is possible 
%
if isnan(out)
    error('Voltage is greater than 10 V or less than 0V. Frequency range is -27.5 to 105.5 MHz');
end

end