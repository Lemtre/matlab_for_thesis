function out = ldpc_create(codename,chkbits,totbits,ldpcmethod,genmethod,seed)

%% set defaults
if nargin < 4 || isempty(ldpcmethod)
  ldpcmethod = 'evenboth 3 no4cycle';
end
if nargin < 5 || isempty(genmethod)
  genmethod = 'dense';
end
if nargin < 6 || isempty(seed)
  seed = 1;
end

%% figure out path of executables
pathstr = fileparts(mfilename('fullpath'));
exe1 = [pathstr filesep 'make-ldpc'];
exe2 = [pathstr filesep 'make-gen'];
randfile1 = [pathstr filesep 'randfile'];
randfile2 = [pwd filesep 'randfile'];

%% generate the LDPC code
if ~exist(randfile2,'file')
  copyfile(randfile1,randfile2);
end
[s,r] = system(sprintf('"%s" "%s.pck" %i %i %i %s',exe1,codename,chkbits,totbits,seed,ldpcmethod));
if s ~= 0
  error(r);
end
out{1} = r;
[s,r] = system(sprintf('"%s" "%s.pck" "%s.gen" %s',exe2,codename,codename,genmethod));
if s ~= 0
  error(r);
end
out{2} = r;
