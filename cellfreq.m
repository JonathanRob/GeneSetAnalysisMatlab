function [item,freq] = cellfreq(c)
%cellfreq  Calculate the frequency of each item in a cell array.
%
% Usage:
%
%   [item,freq] = cellfreq(c);
%
% Input: 
%
%   c      A cell array (vector or matrix).
%
% Output:
%
%   item   A unique list of the items (entries) in cell array C.
%
%   freq   Frequency at which each unique ITEM occurs in cell array C.
%
%
% Jonathan Robinson, 2019-04-14


item = unique(c,'stable');  % determine unique items in C
[~,c] = ismember(c,item);  % convert C to numeric values to increase speed
edges = 0.5:1:(length(item)+0.5);  % define edges between each value of C
freq = histcounts(c,edges)';  % count occurrence of each value of C


