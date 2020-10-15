function [s, history] = sdp_lasso(A, p, nIter, Initialization)
% SDP_LASSO solve the L-p norm optimization problem: 
%       minimize_s ||s'*A*s||_p,  s.t.  |s|==1
%   via semidefenite programming (SDP) and rank-1 regularization. 
% 
% See 'Convex Optimization and Euclidean Geometry' by Dattorro for more details.
% 
% Stop condition: when the relative change in the objective function 
%                 falls below 'tolerance'.
% Maximum iteration number: 'nIter'.
% 'p': which normal to optimize.
% 'Initialization': starting point of the algorithm.
% history = [total cost, original objective, projection, eigenvalue residual]
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

%% Tunable parameters
tolerance = 1e-4;
if length(A) < 10
    zeta = 2;
else
    zeta = 16;
end
%% Initialization
history = NaN(nIter+1, 4);

if nargin<4
    Initialization = 'gauss';
end
dim = length(A{1});
if strcmp(Initialization, 'global')
    [X0, history(1,1)] = sdp_rglrz(A, p, zeros(dim), 0);
elseif strcmp(Initialization, 'gauss')
    X0 = zeros(dim);
    X0(1) = 1;
else
    error('Invalid initialization.');
end
[W, ~, D0] = gen_projmtx_rank1(X0);
z0 = cellfun(@(A_) trace(A_*X0), A);
[history(1,1),history(1,2)] = deal(norm(z0(2:end), p));
history(1,3) = trace(W*X0);
history(1,4) = sum(D0(2:end));

%% Iterations
clear sdp_rglrz
t0 = tic;
figure;
for ii = 1 : nIter
    t1 = tic;
    [X, history(ii+1,1)] = sdp_rglrz(A, p, W, zeta);
    if isempty(X)
        display('Empty output');
        break
    end
    history(ii+1,3) = real(trace(W*X));
    history(ii+1,2) = history(ii+1,1) - zeta*history(ii+1,3); 
    
    %% Generate projection matrix "W" for the next iteration
    [W, V_sorted, D_sorted] = gen_projmtx_rank1(X);
    s = V_sorted(:,1);
    history(ii+1,4) = sum(D_sorted(2:end));
    
    %% Print results and time consumption
    tt0 = toc(t0);
    tt1 = toc(t1);
    fprintf(1, 'Iteration / total time: %0.1f / %0.1f seconds\n', tt1, tt0);
    fprintf(1, 'Original objective: %f, Eigenvalue residual: %f\n',...
        history(ii+1,2),history(ii+1,4));
    %% Plot results
    drawnow;
    plot(1:nIter+1, history, 'LineWidth', 1); axis tight
    legend('Total cost', 'Original objective', 'Projection', 'Eigenvalue residual');    
    %% Check convergency
    if abs(history(ii,1)-history(ii+1,1))<tolerance
        display(['Optimization has converged and terminated at the '...
            num2str(ii) 'th iteration.']);
        break
    end
end
clear sdp_rglrz

%% Generate rank-1 projection matrix
    function [W, V_sorted, D_sorted] = gen_projmtx_rank1(X)
        [V, D] = eig(X, 'vector');
        [D_sorted, idx] = sort(D, 'descend');
        V_sorted = V(:,idx);
        U = V_sorted(:,2:end);
        W = U*U';
    end
end

