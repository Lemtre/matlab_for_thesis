function [a, cost_iter] = iter_opt(A, weights, obj, nIter)
% ITER_OPT iteratively solve the optimization problem: 
%       minimize_a  weights.*(|a'*A*a|-obj).^2.
% See DOI: 10.1109/JOE.2013.2291139 for more details.
% 
% Stop condition: when the relative change in the objective function 
%                 falls below 'tolerance'.
% Maximum iteration number: 'nIter'.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.
A = A(:);
dim = length(A{1});
a = [1; zeros(dim-1,1)];
tolerance = 1e-3;
cost_iter = NaN(nIter+1, 1);
cost_iter(1) = abs(cellfun(@(X_) a'*X_*a, A)).^2' * weights;
gamma = diag(weights);
for ii = 1 : nIter
    B = cell2mat(cellfun(@(X_) a'*X_, A, 'UniformOutput',false));
    a_n = (B'*gamma*B)^-1*B'*gamma*obj;
    a = real(a_n+a)/2;
    err = B*a - obj;
    cost_iter(ii+1) = err'*gamma*err;
    if abs(cost_iter(ii+1)-cost_iter(ii))<tolerance
        display(['Optimization has converged and terminated at the ' num2str(ii) 'th iteration.']);
        return
    end
end

end