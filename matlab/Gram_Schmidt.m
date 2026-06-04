function [Q, R] = Gram_Schmidt(A)
% Gram_Schmidt: QR factorization using Gram-Schmidt (modified style)
% Input:  A is n x p
% Output: Q is n x p with orthonormal columns, R is p x p upper triangular

[n, p] = size(A);

Q = zeros(n, p);
R = zeros(p, p);

for j = 1:p
    x = A(:, j);                 % x is n x 1

    if j == 1
        r = [];                  % no previous columns
        v = x;
    else
        r = Q(:, 1:j-1)' * x;     % r is (j-1) x 1
        v = x - Q(:, 1:j-1) * r;  % v is n x 1
        R(1:j-1, j) = r;
    end

    R(j, j) = norm(v);

    if R(j, j) < 1e-12
        error('Column %d is linearly dependent on previous columns.', j);
    end

    Q(:, j) = v / R(j, j);
end
end