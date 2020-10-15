function Q = include_region_OQAM(K, M, ratio_null, method, varargin)
% INCLUDE_REGION_OQAM find those grid points that need to be constrained
% for OQAM modulation.
% 
%   K:            number of subcarriers.
%   M:            number of subsymbols.   
%   ratio_null:   ratio of nulls-broadening.
%   method:       OQAM modulation type, 'cmt' or 'smt'.
%   varargin:     'edgesOnly', logical variable,
%                 whether to include edges only.
%   Q:            cell array, constrained coordinates.
% 
% Author: Jing Tian
% Last modified: 9th Oct, 2019.
if nargin == 4
    edgesOnly = true;
elseif nargin > 6
    error('The number of optional parameters should be no greater than two.');
elseif nargin == 5
    error('Optional parameters should always go by pairs.');
elseif ~strcmp(varargin{1}, 'edgesOnly')
    prop = varargin{1};
    error(['There is no ' prop ' property of this function.']);
else
    edgesOnly = varargin{2};
end

switch lower(method)
    case 'cmt'
        grids_idx = {[0 1], [2 0], [2 1], [0 2], [4 0]};
    case 'smt'
        grids_idx = {[0 2], [1 0], [1 2], [2 0]};
    otherwise
        error('Non-specific method');
end

xspan = ratio_null/M;
yspan = ratio_null/K;
Q = cellfun(@(idx) get_coordinates(idx, edgesOnly), grids_idx(:),...
    'Uniformoutput', false);
Q = num2cell(cell2mat([[0 0]; Q]), 2);

    function c = get_coordinates(idx, edge_flag)
        ovr = 10; % oversampling rate
        xcent = idx(1)/M;
        ycent = idx(2)/K;
        if ~edge_flag
            xx = unique(linspace(xcent-xspan/2, xcent+xspan/2, ovr));
            yy = unique(linspace(ycent-yspan/2, ycent+yspan/2, ovr));
        else
            xx = unique([xcent-xspan/2, xcent+xspan/2]);
            yy = unique([ycent-yspan/2, ycent+yspan/2]);
        end
        [X,Y] = meshgrid(xx,yy);
        c = [X(:) Y(:)];
    end

end