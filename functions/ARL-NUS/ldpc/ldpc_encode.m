function [ebits,out] = ldpc_encode(codename,bits)

%% figure out path of executables
pathstr = fileparts(mfilename('fullpath'));
exe1 = [pathstr filesep 'encode'];

%% write input data to file
tfn = tempname;
f = fopen(tfn,'wt');
for j = 1:size(bits,2)
  fprintf(f,'%i\n',bits(:,j));
end
fclose(f);

%% run encode on the file
[s,r] = system(sprintf('"%s" "%s.pck" "%s.gen" "%s" "%s-out"',exe1,codename,codename,tfn,tfn));
if s ~= 0
  error(r);
end
out{1} = r;

%% read results and delete intermediate file
f = fopen([tfn '-out'],'rt');
j = 1;
s = fgetl(f);
while ischar(s)
  ebits(:,j) = double(s(:))-48;
  j = j + 1;
  s = fgetl(f);
end
fclose(f);
delete(tfn);
delete([tfn '-out']);
