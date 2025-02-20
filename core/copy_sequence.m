function [snew,opt] = copy_sequence(function_handle,varargin)

function_name = func2str(function_handle);
file_name = which([function_name,'.m']);
s = fileread(file_name);

%Find first line break. This should work on all modern operating systems
r = regexp(s,'\n');
%Break file into to parts, one before line break and one after
s1 = s(1:r(1));
s2 = s((r(1) + 1):end);
%Create new string to insert
sinsert = ['%% These were the input arguments',sprintf('\r\n')];
% if numel(varargin) == 1 && isa(varargin{1},'SequenceOptions')
%     sinsert = [sinsert,sprintf('varargin{1} = %s\r\n',varargin{1}.print)];
% else
    for nn = 1:numel(varargin)
        if isa(varargin{nn},'SequenceOptions')
            sinsert = [sinsert,sprintf('varargin{%d} = %s\r\n',varargin{1}.print)];
        elseif ischar(varargin{nn}) || isstring(varargin{nn})
            sinsert = [sinsert,sprintf('varargin{%d} = %s;\r\n',nn,varargin{nn})];
        else
            sinsert = [sinsert,sprintf('varargin{%d} = %.6g;\r\n',nn,varargin{nn})]; %#ok<*AGROW>
        end
    end
% end
%Insert string into file
snew = [s1,sinsert,s2];

opt = parse_maker_variable_argument_list(varargin{:});

