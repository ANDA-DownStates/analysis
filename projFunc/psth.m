function [ output ] = psth( spikeTrain, range, binSize )
%return post stimulation time histogram
%   range is a 2 elements vector [range_start, range_end]
%   This function will reject spikes between the last edge to the right
%       range in case the gap between the range is not a multiple of bin.

edges = range(1):binSize:range(2);
% edges = [edges, range(2)];    % exclude spikes between the last edge and the real right range.
output = histcounts(spikeTrain, edges);

end

