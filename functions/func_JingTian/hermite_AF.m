function [x_out, y_out, AF] = hermite_AF(syn_para1, syn_para2,...
    x_range, y_range, coordinate_type)
% HERMITE_AF calculate the closed-form ambiguity function between two
%   signals synthesized with Hermite functions (HFs).
%
%   syn_para1:       the orders and coefficients of HFs to synthesize signal #1
%   syn_para2:       the orders and coefficients of HFs to synthesize signal #2
%   x_range:         the range of angle in polar coordinates/the range of delay in Cartesian coordinates
%   y_range:         the range of radius in polar coordinates/the range of frequency in Cartesian coordinates
%   coordinate_type: in Cartesian coordinates or polar coordinates
%   [x_out, y_out]:  the output mesh grids in Cartesian coordinates
% 
% Example:
%   coordinate_type = 'polar';
%   x_range = linspace(-pi, pi, 180);
%   y_range = linspace(0, 6);
%   ord1 = [0; 4]; ord2 = [0; 4];
%   coe1 = [1; 0.7]; coe2 = [1; 0.2];
%   syn_para1 = [ord1 coe1];
%   syn_para2 = [ord2 coe2];
%   [x,y,z] = hermite_AF(syn_para1, syn_para2, x_range, y_range, coordinate_type);
%   z = abs(z); z = z./max(max(z)); 
%   figure; contour(x,y,abs(z),100)
% 
% 
% This function uses the closed-form expression from
% DOI:10.1109/TSP.2015.2488580,
% the calculation might take some time since it involves symbolic
% functions.
% 
% Author: Jing Tian
% Last modified: 9th Sep, 2020.

if strcmp(coordinate_type, 'polar')
    [tt, rr] = meshgrid(x_range, y_range);
    [x_out, y_out] = pol2cart(tt, rr);
elseif strcmp(coordinate_type, 'cartesian')
    [x_out, y_out] = meshgrid(x_range, y_range);
    [tt, rr] = cart2pol(x_out, y_out);
else 
    error('Unsupported coordinate system.');
end
theta = tt(:);
rr = rr(:);

[order1, coe1, order2, coe2] = deal(syn_para1(:,1), syn_para1(:,2),...
    syn_para2(:,1), syn_para2(:,2));
if sum(fix(order1)~=order1) || sum(fix(order2)~=order2)
    error('Orders of the Hermite functions should be integers.');
end
ords = meshQ(order1, order2);
coes = meshQ(coe1, coe2);
ord_sum = sum(ords, 2);
ord_cell = num2cell(ords, 2).';
coe_cell = num2cell(coes, 2).';
[ord_sum_uniq, ~, ic] = unique(ord_sum);
id_uniq_even = find(mod(ord_sum_uniq,2)==0);
id_uniq_odd = find(mod(ord_sum_uniq,2)~=0);

p_even = arrayfun(@(id) ~isempty(intersect(id, id_uniq_even)), ic);
p_odd = arrayfun(@(id) ~isempty(intersect(id, id_uniq_odd)), ic);
id_n_m = 1:length(ord_sum);
id_even = id_n_m(p_even==1);
id_odd = id_n_m(p_odd==1);

H_ord = NaN(length(rr), length(ord_sum));
if ~isempty(id_even)
    sum_max_even = ord_sum_uniq(id_uniq_even(end));
    k_max_even = floor(ord_sum_uniq(id_uniq_even(end))/2);
    kk_even = 0 : k_max_even;
    h_table_even = arrayfun(@(kk) hermite_func(sum_max_even-2*kk, rr/sqrt(2)),...
        kk_even, 'Uniformoutput', false);
    h_table_even = cell2mat(h_table_even);
    H_ord_even = cellfun(@(ord) sum_ovr_h(ord(1), ord(2), h_table_even),...
        ord_cell(id_even), 'Uniformoutput', false);
    H_ord(:, id_even) = cell2mat(H_ord_even);
end

if ~isempty(id_odd)
    sum_max_odd = ord_sum_uniq(id_uniq_odd(end));
    k_max_odd = floor(ord_sum_uniq(id_uniq_odd(end))/2);
    kk_odd = 0 : k_max_odd;
    h_table_odd = arrayfun(@(kk) hermite_func(sum_max_odd-2*kk, rr/sqrt(2)),...
        kk_odd, 'Uniformoutput', false);
    h_table_odd = cell2mat(h_table_odd);
    H_ord_odd = cellfun(@(ord) sum_ovr_h(ord(1), ord(2), h_table_odd),...
        ord_cell(id_odd), 'Uniformoutput', false);
    H_ord(:, id_odd) = cell2mat(H_ord_odd);
end
H_cell = num2cell(H_ord, 1);
H_syn = cellfun(@(coe, ord, H) coe(1)*conj(coe(2))*(-1j)^sum(ord)*sqrt(pi).*exp(1j*theta*diff(ord)).*H,...
    coe_cell, ord_cell, H_cell, 'Uniformoutput', false);
H_sum_ovr = sum(cell2mat(H_syn), 2);
AF = reshape(H_sum_ovr, [], length(x_range));


function c = coe_calc(k,n,m)

uceil = min([k m n]);
u_exprs = @(u)(-4).^u .* factorial(n+m-2*u)./(factorial(u).*factorial(m-u)...
    .*factorial(n-u).*factorial(k-u));
u = sum(u_exprs(0 : uceil), 2);
c = sqrt(factorial(n)*factorial(m))*u/(2.^k*sqrt(2^(n+m)*factorial(n+m-2*k)*sqrt(pi)));

end

function grids = meshQ(x ,y)
    [X, Y] = meshgrid(x, y);
    grids = [X(:) Y(:)];
end

function H = sum_ovr_h(n, m, h_table)
kk = 0 : floor((m+n)/2);
c_k = arrayfun(@(k) coe_calc(k,n,m), kk);
H = repmat(c_k, size(h_table,1), 1).* h_table(:, end-kk(end)+kk);
H = sum(H,2);
end

end