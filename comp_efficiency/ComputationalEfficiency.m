function [efficiency, gammabar] = ComputationalEfficiency(p, m, n)

if nargin < 2
    m = 100;
end
if nargin < 3
    n = 30*m;
end

gammabar = (1:m)/m;
b = NumBatches(gammabar,m,n);

% If number of batches doesn't change efficiency won't change either.
% Remove duplicate batch numbers
% for ii=length(gammabar):-1:2
%     if b(ii-1) == b(ii)
%         gammabar(ii-1) = [];
%         b(ii-1) = [];
%     end
% end

T = gammabar.^p;

efficiency = 1./(T.*NumBatches(gammabar,m,n).*RelativeVariance(gammabar));
efficiency = efficiency/efficiency(end);