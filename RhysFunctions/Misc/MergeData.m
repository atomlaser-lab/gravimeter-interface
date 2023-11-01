function [MergedData] = MergeData(PreviousData,NewData)

% % This works when scanning 1 variable

if strcmpi(PreviousData.ParamName,NewData.ParamName) ~= 1
    warning('Merged Data may be different quantities')
end
MergedData.ParamName = PreviousData.ParamName;
MergedData.ParamUnits = PreviousData.ParamUnits;

MergedData.param = [PreviousData.param,NewData.param];

MergedData.N = [PreviousData.N;NewData.N];
MergedData.peakOD = [PreviousData.peakOD;NewData.peakOD];
MergedData.Nsum = [PreviousData.Nsum;NewData.Nsum];
MergedData.files = [PreviousData.files;NewData.files];
MergedData.becFrac = [PreviousData.becFrac;NewData.becFrac];
MergedData.T = [PreviousData.T;NewData.T];
MergedData.R = [PreviousData.R;NewData.R];
MergedData.xPos = [PreviousData.xPos;NewData.xPos];
MergedData.yPos = [PreviousData.yPos;NewData.xPos];
MergedData.xwidth = [PreviousData.xwidth;NewData.xwidth];
MergedData.ywidth = [PreviousData.ywidth;NewData.ywidth];
MergedData.BECWidthX = [PreviousData.BECWidthX;NewData.BECWidthX];
MergedData.BECWidthY = [PreviousData.BECWidthY;NewData.BECWidthY];
end