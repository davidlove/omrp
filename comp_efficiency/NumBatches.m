function n = NumBatches( gammabar, m, n )
%NUMBATCHES computes number of batches for gammabar

if nargin < 2
    m = 100;
end

if nargin < 3
    n = 30*m;
end

n = floor((n-m)./(m.*gammabar)) + 1;

end
