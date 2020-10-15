%% 
clear
close all
dbstop if error
%% Set parameters
set_para_HGLfiltsyn;
%% Generate dilated HGL pairs
g0 = gen_Pei_HGLvec0(N_eig, floor(N_eig/2)/2); % generate 0th order HGL vector
H_eig = gen_HGLeigvec(g0, orders); % generate 4n-th order HGL vectors
H = per_upsamp(H_eig, N_gfdm); % upsample to get scaled HGL vectors 
H = fftshift(H,1); % shift to center
%% Compute Aq matrices 
Q = include_region_OQAM(K, M, ratio_null, method);
% ----------------------------------------------------------------------
A = compute_A_mtx(H, Q); % compute A matrices
% ----------------------------------------------------------------------
%% Solution #1: Iterative optimization
obj = [1; zeros(length(A)-1,1)]; % objective
weights = [1; 40*ones(length(A)-1, 1)]; % weights
% ----------------------------------------------------------------------
% iterative optimization
[a1, cost_iter] = iter_opt(A, weights, obj, nIter);
%% Solution #2: SDP relaxation + LASSO 
% only including the edges will greatly speed up the LASSO algorithm
Q_lasso = include_region_OQAM(K, M, ratio_null, method, 'edgesOnly', true);
A_lasso = compute_A_mtx(H, Q_lasso); 
[a2, history] = sdp_lasso(A, p, nIter, Initialization);

%% Solution #3: Real mismatched filter pair design
% different Tx and Rx filters
[a3, cost_asym] = asymm_real_altopt(A, nIter);
%% synthesize filter
g1 = H*a1; 
g1s = ifftshift(g1, 1); % prototype filter(s)

g2 = H*a2;
g2s = ifftshift(g2, 1);

g3 = H*a3;
g3s = ifftshift(g3, 1);
%% Plot results
display('Plotting results...');
% ----------------------------------------------------------------------
% Solution #1:
[paf1, At1, Af1] = plot_result_filtdes(g1, K, M, ratio_null, method);
check_OQAMmtx(g1s, K, M, method);
% ----------------------------------------------------------------------
% Solution #2:
[paf2, At2, Af2] = plot_result_filtdes(g2, K, M, ratio_null, method);
check_OQAMmtx(g2s, K, M, method);
% ----------------------------------------------------------------------
% Solution #3:
[paf3, At3, Af3] = plot_result_filtdes(g3, K, M, ratio_null, method);






