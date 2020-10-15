function [symint_demod, bit_out] = bbSymDemap(modType, sym)
% BBSYMDEMAP map complex symbols to integers and 0-1 bits via hard decision.
% 
%   modType:      modulation type.
%   sym:          symbols.
%   symint_demod: integers after remapping.
%   bit_out:      0-1 bits after demodulation.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

if ismatrix(sym) && ~iscell(sym)
    sym_cell = num2cell(sym, 1);
end

switch modType
    case 'BPSK'
        modOrder = 2;
        demapFunc = @pskdemod;
    case 'QPSK'
        modOrder = 4;
        demapFunc = @pskdemod;
    case '8PSK'
        modOrder = 8;
        demapFunc = @pskdemod;
    case '16QAM'
        modOrder = 16;
        demapFunc = @qamdemod;
    otherwise
        warning('Invalid symbol mapping scheme.');
        return;
end
ini_phase = 0;
symOrder = 'gray';

symint_cell = cellfun(@(sym_col) demapFunc(sym_col, modOrder, ini_phase, symOrder),...
    sym_cell, 'Uniformoutput', false);
bit_cell = cellfun(@(x) de2bi(x), symint_cell, 'Uniformoutput', false); % convert decimal numbers to binary bits
bit_out = cell2mat(bit_cell); % convert cell array to matrix
bit_out = reshape(bit_out, size(sym,1), [], size(sym,2));
symint_demod = cell2mat(symint_cell);


end