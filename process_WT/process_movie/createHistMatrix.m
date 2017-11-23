function [ x_bin_centers, y_bin_centers,counts ] = createHistMatrix( scores1, scores2,nBins )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nBins_x = linspace(min(scores1),max(scores1),nBins);
nBins_y = linspace(min(scores2),max(scores2),nBins);

[counts, bin_centers] = hist3([scores1,scores2], {nBins_x nBins_y});
x_bin_centers = bin_centers{1};
y_bin_centers = bin_centers{2};
end

