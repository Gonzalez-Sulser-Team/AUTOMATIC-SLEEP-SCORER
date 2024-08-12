function [ind]=find_in_interval(intervals, ts)
% [ind]=find_in_interval(intervals, timestamps)
% function to find parts of a vector (eg event timestamps) that fall into
% certain intervals (eg sleep epochs).
% inputs :
%     intervals - matrix of pairs defining edges (start/stop) of intervals
%                 eg:    intervals = [ 1 10; 20 25];
%     process   - vector of timestamps or else
%                 eg process=[ 0 3 5 15 16 21]
% outputs:
%     ind       - logical index vector of size process 
%                 eg [0 1 1 0 0 1]=find_in_interval(intervals, process)
% units of time must be conistent

% (c) Ullrich Bartsch, University of Bristol 2014
% ubartsch(at)gmail.com
%% --------------------------------------------------------------------------

ind=false(size(ts));
for i=1:size(intervals,1)
    ind( ts>intervals(i,1) & ts<intervals(i,2) )= true;
end


