function [H_BP, id] = syn_chanmtx_BP(xi, chanmtx)
% SYN_CHANNELMTX synthesize the channel matrix with BP estimation result.
% 
%   xi:      solution of the BP optimization problem.
%   chanmtx: cell array containing channel matrices.
%   H:       frequency-domain channel matrix.
% 
% Authour: Jing Tian
% Last modified: 10th Sept, 2020.

xi = xi(:);
rho = 0.01;
xi_max = max(abs(xi));
id = find(abs(xi)>=rho*xi_max); % discard entries whose energy is -20dB lower
xi_none0 = num2cell(xi(id), 2);
Hps = chanmtx(id);

H_BP = cellfun(@(a,H) a*H, xi_none0, Hps(:), 'Uniformoutput', false);
H_BP = full(cell2mat(H_BP(:).'));
H_BP = reshape(H_BP, [size(Hps{1}) length(id)]);
H_BP = sum(H_BP,3);

end