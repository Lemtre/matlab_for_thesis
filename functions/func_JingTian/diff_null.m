function [diff_sym, idx1] = diff_null(symvec, null_pos)
% DIFF_NULL differentially demodulate the input symbol vector, exclude
% all the null carriers and those singleton carriers between two nulls.
% 
%   symvec:   OFDM symbols after FFT, vector.
%   null_pos: position of null carriers, vector.
%   diff_sym: differentially demodulated symbols.
%   idx1:     positions of the diff-demodulated symbols in the original vector.
% 
% Input SYM should be vector instead of matrix, because null subcarriers
% are randomly positioned for different OFDM blocks, therefore there could
% be different numbers of differentially decoded symbols, causing dimensions 
% mismatch for matrix operation.
% 
% Author: Jing Tian
% Last Modified: 25 May 2018

    if numel(symvec) > length(symvec) || numel(null_pos) > length(null_pos)
        error('The first input should be a vector!');
    elseif length(null_pos) > length(null_pos)
        error('There are more null carriers than symbols!');   
    end
    symvec = symvec(:);
    null_pos = null_pos(:);
    
    filter_all = true(size(symvec));
    singleton_pos = null_pos(diff([0; null_pos]) == 2) - 1; % Find positions of singleton carriers between two nulls.
    pos = union(null_pos, singleton_pos);
    filter_all(pos) = false;
    filter_ini = filter_all(1:end-1);
    filter_end = filter_all(2:end);

    find_interval = diff(filter_all);

    filter_ini(find_interval == -1) = false;
    filter_end(find_interval == 1) = false;

    len = (1:length(symvec)).';
    idx1 = len(filter_ini);
    idx2 = len(filter_end) + 1;

%     diff_sym = sym(idx1) .* conj(sym(idx2));
    diff_sym = symvec(idx1) ./ symvec(idx2);
    
    if size(symvec,1) < size(symvec,2) % If the input is a row vector, so is the output
        diff_sym = diff_sym.';
    end
    
end