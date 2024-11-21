function makeFunctionSignatures

sq = initSequence;

N = sq.numChannels;
names = cell(N,1);
for nn = 1:N
    names{nn} = sq.channels(nn).name;
end

choices = '';
for nn = 1:numel(names)
    choices = sprintf('%s''%s'',',choices,names{nn});
end
choices = sprintf('choices={%s}',choices(1:end-1));

s = fileread('core/functionSignatures.json');

r = regexprep(s,'choices\=\{.*?\}',choices,'dotexceptnewline');

fid = fopen('core/functionSignatures.json','w');
fprintf(fid,r);
fclose(fid);


end