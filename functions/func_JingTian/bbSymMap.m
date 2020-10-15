function [sym_cplx, sym_int] = bbSymMap(bit_in, modType)
% BBSYMMAP map 0-1 bits to complex symbols.
% 
%   bit_in:   input bit array,
%             1-dim - number of symbols per block
%             2-dim - modulation order
%             3-dim - number of blocks
%   sym_cplx: complex symbols drawn from a certain constellation
%   sym_int:  integers drawn from a certain alphabet
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.


if ~ismatrix(bit_in) || numel(size(bit_in))>3
    error('The first input should be an array of dimensions less or equal than three.');
end

switch modType
    case 'BPSK'
        modOrder = 2;
        mapFunc = @pskmod;
    case 'QPSK'
        modOrder = 4;
        mapFunc = @pskmod;
    case '8PSK'
        modOrder = 8;
        mapFunc = @pskmod;
    case '16QAM'
        modOrder = 16;
        mapFunc = @qammod;
    otherwise
        warning('Invalid symbol mapping scheme.');
        return;
end
ini_phase = 0;
symOrder = 'gray';

if log2(modOrder)~=size(bit_in,2)
    error('Data dimension mismatch!');
end
bit_cell = squeeze(num2cell(bit_in, [1 2])).'; % convert array to cell array
sym_int_cell = cellfun(@(x) bi2de(x), bit_cell, 'Uniformoutput', false); % convert binary data to decimal numbers
sym_int = cell2mat(sym_int_cell); % convert cell array to matrix

sym_cplx = mapFunc(sym_int, modOrder, ini_phase, symOrder);

end