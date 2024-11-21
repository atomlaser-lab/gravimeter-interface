function [VoltageNeeded] = RFtoV(DesiredFrequency)

if any(DesiredFrequency > 20 | DesiredFrequency < 0)
    error('Frequency Range is 20 to 0 MHz')
end

% VoltageNeeded = (DesiredFrequency-10)/2;
VoltageNeeded = (DesiredFrequency - 10.0193)/1.9902;

end



