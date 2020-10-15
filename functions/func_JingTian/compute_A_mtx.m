function A = compute_A_mtx(H, Q)

dim = size(H, 2);
[idm, idn] = meshgrid(1:dim, 1:dim);
id_mn = num2cell([idm(:) idn(:)], 2);

display('----------------------------------------------------------------');
display('Computing A matrices...');
t0 = tic;
A = cellfun(@(q) compute_Aq(q), Q, 'Uniformoutput', false);
t1 = toc(t0);
display('A matrices are ready!');
fprintf(1, 'Total elapsed time: %0.1f seconds\n', t1);
display('----------------------------------------------------------------');

    function Aq = compute_Aq(q)
        Aqs = cellfun(@(id) paf_eqv(H(:,id(1)),H(:,id(2)),q(1),q(2)),id_mn);
        Aq = reshape(Aqs, dim, dim);
    end

end