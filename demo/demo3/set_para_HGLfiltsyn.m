%% Set GFDM/OQAM parameters
method = 'smt'; % OQAM method
% method = 'cmt';
% ----------------------------------------------------------------------
rho = 2; % density
% ----------------------------------------------------------------------
M = 4; % number of subsymbols
K = 32; % number of subcarriers
if strcmp(method, 'cmt') && M<8
    warning('Block length is too short for CMT method. Please increase M.');
end

switch lower(method)
    case 'cmt'
        N_eig = M^2/rho;
    case 'smt'
        N_eig = M^2*rho;
    otherwise
        error('Non-specified method');
end
N_gfdm = K*M;
% clearvars -except K M N_eig N_gfdm;
%% Number of HGL vectors to include
dim = 4;
orders = 4*(0:dim-1); % 4n-th order HGL vectors
%% nulls broadening region parameters
ratio_null = 0.1;
%% Maximum number of iterations
% % Iterative algorithms usually converge within tens of iterations
nIter = 100;
%% which norm to optimize
p = 2;
% p = inf;
%% rank-1 LASSO initialization
% choose starting point:
% Initialization = 'global'; % global optimizer of SDP relaxation
Initialization = 'gauss'; % 0th HGL vector, more robust.  






