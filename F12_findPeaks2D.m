function loc = findPeaks2D(im,quality)
% Returns local maxima in the image. quality is a threshold
% expressed as a fraction of the maximum value of metric.

%   Copyright 2010-2014 The MathWorks, Inc.

maxMetric = max(im(:));
if maxMetric <= eps(0)
    loc = zeros(0,2); %empty matrix
else
    bw = imregionalmax(im, 8);
    threshold = quality * maxMetric;
    bw(im < threshold) = 0;
    bw = bwmorph(bw, 'shrink', Inf);
    
    % Exclude points on the border
    bw = imclearborder(bw,8);
    
    % Find location of the peaks    
    idx = find(bw);
    loc = zeros([length(idx) 2]);
    [loc(:, 2), loc(:, 1)] = ind2sub(size(im), idx);
end