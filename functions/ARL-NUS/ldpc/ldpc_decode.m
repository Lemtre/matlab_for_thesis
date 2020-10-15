function [bits,out] = ldpc_decode(codename,ebits,ch,method)

%% defaults
if ~ischar(ch)
  ch = sprintf('bsc %0.5f',ch);
end
if nargin < 4 || isempty(method)
  method = 'prprp 250';
end
if ~ischar(method)
  method = sprintf('prprp %i',method);
end

%% figure out path of executables
pathstr = fileparts(mfilename('fullpath'));
exe1 = [pathstr filesep 'decode'];
exe2 = [pathstr filesep 'extract'];

%% write input data to file
tfn = tempname;
f = fopen(tfn,'wt');
for j = 1:size(ebits,2)
  fprintf(f,'%i\n',ebits(:,j));
end
fclose(f);

%% run encode on the file
[s,r] = system(sprintf('"%s" "%s.pck" "%s" "%s-out" %s %s',exe1,codename,tfn,tfn,ch,method));
if s ~= 0
  error(r);
end
out{1} = r;
[s,r] = system(sprintf('"%s" "%s.gen" "%s-out" "%s"',exe2,codename,tfn,tfn));
if s ~= 0
  error(r);
end
out{2} = r;

%% read results and delete intermediate file
f = fopen(tfn,'rt');
j = 1;
s = fgetl(f);
while ischar(s)
  bits(:,j) = double(s(:))-48;%#ok
  j = j + 1;
  s = fgetl(f);
end
fclose(f);
delete(tfn);
delete([tfn '-out']);
