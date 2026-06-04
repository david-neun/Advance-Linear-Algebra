function [A, eigApprox, A_hist] = QR_Algorithm(A, numIter)
% QR_Algorithm: unshifted QR iteration using Gram_Schmidt from Problem 7
% Inputs:
%   A        - n x n matrix (should be square)
%   numIter  - number of iterations
% Outputs:
%   A        - final iterate A_numIter
%   eigApprox- diagonal of final A (eigenvalue approximations)
%   A_hist   - (optional) history of iterates

[n, m] = size(A);
if n ~= m
    error('A must be square.');
end

A_hist = zeros(n, n, numIter + 1);
A_hist(:, :, 1) = A;

for k = 1:numIter
    [Q, R] = Gram_Schmidt(A);  % function from Problem 7

    % next similar matrix
    A = R * Q;

    A_hist(:, :, k + 1) = A;
end

eigApprox = diag(A);
end
