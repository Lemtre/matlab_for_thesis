function Y = per_upsamp(X, N)

delay = 1/N: 1/N : 1-1/N;
Xc = num2cell(X,1);
Xd = cellfun(@(h) per_interp_delay(h,delay), Xc, 'Uniformoutput', false);
Y = cellfun(@(hd,h0) reshape([h0(1) hd(1,:)].', [], 1), Xd, Xc, 'Uniformoutput', false);
Y = cell2mat(Y);

end