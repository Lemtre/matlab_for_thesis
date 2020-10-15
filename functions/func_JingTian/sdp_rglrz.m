function [X_opt, optval] = sdp_rglrz(A, p, W, zeta)

dim = length(A{1});

% cvx_solver mosek
cvx_begin sdp
    cvx_precision high;
    variable X(dim,dim) symmetric;
    expression z(length(A),1);
    for ii = 1:length(A)    
        z(ii) = trace(A{ii} * X);%#ok
    end
    minimize(norm(z(2:end), p) + zeta*trace(W*X));
    subject to
        X == semidefinite(dim);%#ok
        z(1) == 1;%#ok
cvx_end

X_opt = X;
optval = cvx_optval;

end