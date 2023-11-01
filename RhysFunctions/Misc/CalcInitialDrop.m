function DropTime = CalcInitialDrop(DesiredSize,InitialSize,Temp)
% Altin Eq: 3.9
% sig(t)^2 = sig(0)^2 + t^2*sig_v^2
% sig_v = sqrt(2*kb*T/m)

% This function uses Equantion 3.9 seen in Altin's thesis to calculate the
% expansion time required to obtain a cloud with a desired size provided
% the initial size and temperature is known.

sig_v = sqrt(2*const.kb*Temp/const.mRb);
DropTime = sqrt((DesiredSize^2 - InitialSize^2)/(sig_v^2));


end