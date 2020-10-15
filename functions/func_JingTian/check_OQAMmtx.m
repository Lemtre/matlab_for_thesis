function check_OQAMmtx(g, K, M, method)
% CHECK_OQAMMTX check the real orthogonality between OQAM matrices.
g = g(:);
if length(g) ~= K*M
    error('Filter length mismatch.');
end
[Ai, Aq] = gen_GFDM_OQAMmtx(g, K, M, 'sc', method);
figure; set(gcf,'unit','centimeters','position',[27 42 30 30]);
subplot(2,2,1); imagesc(real(Ai'*Ai)); title('\Re\{Ai^HAi\}'); axis square;
subplot(2,2,2); imagesc(real(Aq'*Aq)); title('\Re\{Aq^HAq\}'); axis square;
subplot(2,2,3); mesh(real(Ai'*Aq)); axis tight; title('\Im\{Ai^HAq\}');
subplot(2,2,4); mesh(real(Aq'*Ai)); axis tight; title('\Im\{Aq^HAi\}');

end