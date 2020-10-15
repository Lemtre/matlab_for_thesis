function bit = gen_binbits(modType, Nsym, Nblk)
% GEN_BINBITS generate 0-1 binary bit stream.
% 
%   modType: symbol mapping scheme.
%   Nsym:    number of symbols per block.
%   Nblk:    number of blocks.
%   bit:     output bits.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.


switch modType
    case 'BPSK'
        modOrder = 2;
    case 'QPSK'
        modOrder = 4;
    case '8PSK'
        modOrder = 8;
    case '16QAM'
        modOrder = 16;
    otherwise
        warning('Invalid symbol mapping scheme.');
        return;
end
bit = randi([0 1], Nsym, log2(modOrder), Nblk); % generate 0-1 bit stream

end