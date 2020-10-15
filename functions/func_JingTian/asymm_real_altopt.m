function [a, cost] = asymm_real_altopt(A, nIter)
% ASYMM_REAL_ALTOPT solve the optimization problem, 
%       minimize_{x,y}  |x'*A*y|^2,
%  via alternate minimization. The key idea is similar to 
%       DOI: 10.1109/ICASSP.2011.5946241,
%  other than constraining 'x' and 'y' to be real in this function. 
% 
% Stop condition: when the relative change in the objective function 
%                 falls below 'tolerance'.
% Maximum iteration number: 'nIter'.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

A = A(:);
% A1 = real(A{1});
A1 = A{1};
A0 = A(2:end);
%% Initialization
tolerance = 1e-3;

dim = length(A1);
x = [1; zeros(dim-1,1)];
y = [1; zeros(dim-1,1)];
cost = NaN(nIter+1,1);
cost(1) = norm(cellfun(@(AA) y'*AA*x, A0))^2;
%% Alternate minimization
tic
for ii = 1 : nIter
    Rx = update_Rx(y);
    A1_yy = real(A1*(y*y.')*A1');
    x = genEigValProb(Rx, A1_yy);
    
    Ry = update_Ry(x);
    A1_xx = real(A1'*(x*x.')*A1);
    y = genEigValProb(Ry, A1_xx);
    
    a = [x y];
   
% % ----------------------- check convergency ------------------------- % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %            
    cost(ii+1) = y.'*Ry*y;
    if cost(ii+1)-cost(ii)>tolerance
        error(['Objective increases at the ' num2str(ii) 'th iteration.']);
    elseif abs(cost(ii+1)-cost(ii))<tolerance
        display(['Optimization has converged and terminated at the '...
            num2str(ii) 'th iteration.']);
        break
    end
    
end
toc

    function Rx = update_Rx(y)
        Rx = cellfun(@(K_) reshape(K_*(y*y.')*K_', 1, []), A0, 'UniformOutput', false);
        Rx = reshape(sum(real(cell2mat(Rx))), dim, dim);
    end

    function Ry = update_Ry(x)
        Ry = cellfun(@(K_) reshape(K_'*(x*x.')*K_, 1, []), A0, 'UniformOutput', false);
        Ry = reshape(sum(real(cell2mat(Ry))), dim, dim);
    end

    function vmin = genEigValProb(As,Bs)
        epsilon = 1e-5;
%         epsilon = eps;
        As = real(As+As')/2;
        Bs = real(Bs+Bs')/2;
        Bs = Bs + epsilon*eye(dim); % diagonal loading is necessary for stability
        [V,D] = eig(As,Bs);
        [~, ranked] = sort(real(diag(D)),'ascend');
        V_sorted = V(:,ranked);
        vmin = V_sorted(:,1);
    end

end