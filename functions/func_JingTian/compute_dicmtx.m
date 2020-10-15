function A  = compute_dicmtx(ofdmPara, pilots, load_chanmtx, chanmtx_var)
% COMPUTE_DICMTX compute dictionaries for BP channel estimation.
% 
%   ofdmPara:     OFDM parameters.
%                 ofdmPara.idp - pilot positions;
%                 ofdmPara.idobz - positions of observations along the column of A.
%   pilots:       known OFDM pilots.
%   load_chanmtx: logical variable.
%                 true - load channel matrices from chanmtx_var.
%                 false - channel matrices set is already in the workspace.
%                         Pass its value via chanmtx_var.
%   chanmtx_var:  path and filename of the MAT file that contains channel matrices.
%   A:            dictionary matrices.
% 
% Channel matrices should be pre-computed and stored before calling
% COMPUTE_DICMTX.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.


idp = ofdmPara.idp;
idobz = ofdmPara.idobz;

if load_chanmtx
    fprintf(1, '-------------------------------------------------------------\n');
    fprintf(1, 'Loading channel matrices...\n');
    tic;
    var_struct = load(chanmtx_var);
    var_name = char(fieldnames(var_struct));
    Hp = var_struct.(var_name);
    toc;
else
    Hp = chanmtx_var;
end
fprintf(1, 'Computing dictionary matrices...\n');
tic;
A = cellfun(@(H) H(idobz,idp)*pilots(:), Hp, 'Uniformoutput', false);
toc;
A = cell2mat(A(:).');
fprintf(1, 'Dictionary matrices are ready!\n');
fprintf(1, '-------------------------------------------------------------\n');


end