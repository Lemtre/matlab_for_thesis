%% Common settings
modType = 'QPSK';
typeGI = 'cp';
group_order = 'sc';
%% Set bandwidth, center frequeNy, and guard time, etc.
% -------------------------------------------------------------------------
K = 256; % number of subcarriers
M = 4;   % number of subsymbols
N = M*K; % total number of samples per block
Ng = round(N/8); % samples in guard interval
bw = 6e3; % bandwidth
fc = 12e3; % center frequency
fs = 96e3; % sampling frequency of the transmitted signal in passband

sysPara.fc = fc;
sysPara.Tg = Ng/bw;
% -------------------------------------------------------------------------

%% channel parameters
chanPara.Npa = 7; % the number of paths,
chanPara.interarrv_mean = 2e-3; % the mean of inter-arrvals (in seconds)
chanPara.maxAvgPathGainLoss = -20; % average path gain loss (in dB)
chanPara.TimeVar = 1; % channel is time-varying
chanPara.UniDopplerFactor = 0; % path-specific Doppler
chanPara.sigma_v = 0.05; % motion speed (m/s)

snr = 15; 
%% Set positions of pilots
B_coh = 1/sysPara.Tg;
F = bw/N;
lambda = floor(B_coh/F/M);

P = lambda*M;
idp = (1:P:N).';
Np = length(idp);
%% GFDM parameters
var_struct = load('rcos_M4_K256.mat'); % raised cosine filter with roll-off factor 0.1
var_name = char(fieldnames(var_struct));
g_gfdm = var_struct.(var_name);

gfdmPara.filter = g_gfdm;
gfdmPara.Nsc = K;
gfdmPara.Nss = M;
gfdmPara.group_order = group_order;
gfdmPara.id_pilot = idp;

id_AD_sc = 0:lambda:K-1;
id_AD_sc = repmat(id_AD_sc, 1, 2);
id_AD_ss = randi([0 M-1], Np, 1);
id_AD_ss = [id_AD_ss; mod(id_AD_ss+1,M)];
posAD_gfdm = [id_AD_sc(:) id_AD_ss(:)];

%% CFBMC/OQAM parameters
var_struct = load('HGL_M4_K256_smt.mat'); % HGL filter 
var_name = char(fieldnames(var_struct));
g_oqam = var_struct.(var_name);
oqam_method = 'smt';

oqamPara.filter = g_oqam;
oqamPara.Nsc = K;
oqamPara.Nss = M;
oqamPara.group_order = group_order;
oqamPara.oqam_method = oqam_method;
oqamPara.id_pilot = idp;

[id_AD_sc_i, id_AD_sc_q] = deal(0:lambda:K-1);
id_AD_ss_i = 0*ones(Np,1);
id_AD_ss_q = 2*ones(Np,1);
posAD_oqam{1} = [id_AD_sc_i(:) id_AD_ss_i(:)];
posAD_oqam{2} = [id_AD_sc_q(:) id_AD_ss_q(:)];






