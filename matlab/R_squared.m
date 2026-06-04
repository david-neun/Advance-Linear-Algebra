function [R_sq] = R_squared(X, y, beta_prime)
% David Cortes
% March 30, 2026
%
% R_squared - Computes the coefficient of determination, R^2, for a
% least-squares model. R^2 measures how well the model explains the
% variability in the data. A value close to 1 indicates a good fit,
% while a value close to 0 indicates a poor fit.
%
% The function automatically detects whether the ones vector 1 is in
% CS(X) and selects the appropriate SST formula:
%   - If 1 in CS(X): SST = y' * C * y, where C = I - (1/m)*1*1' (centering matrix)
%   - If 1 not in CS(X): SST = y' * y
%
% USAGE:
%   [R_sq] = R_squared(X, y, beta_prime)
%
% INPUTS:
%   X           - Design matrix (m x n)
%   y           - Observed data vector (m x 1 column vector)
%   beta_prime  - Least-squares solution vector (n x 1 column vector)
%
% OUTPUT:
%   R_sq        - Coefficient of determination, R^2 (scalar)
%
% TEST CASES:
%
%   Test 1: Problem 3 - Cricket data (1 in CS(X), polynomial fit)
%     load('cricket_data.mat');
%     chirp_rate = cricket_data(:,2);    % independent variable (m x 1)
%     temperature = cricket_data(:,1);   % dependent variable (m x 1)
%     poly_deg = 1;
%     m = length(chirp_rate);
%     n = poly_deg + 1;
%     X = zeros(m, n);
%     for j = 1:n
%         X(:,j) = chirp_rate.^(j-1);   % Vandermonde matrix (m x n)
%     end
%     beta_prime = X \ temperature;
%     R_sq = R_squared(X, temperature, beta_prime)
%     % Expected: R_sq close to 1 (tests path (h), 1 is in CS(X))
%
%   Test 2: Problem 4 - Metabolic rate data (log-log transform, 1 in CS(X))
%     load('metabolic_rate_data.mat');
%     mass = metabolic_rate_data(:,1);           % independent variable (m x 1)
%     met_rate = metabolic_rate_data(:,2);       % dependent variable (m x 1)
%     X = [ones(length(mass),1), log(mass)];     % log transform on x (m x 2)
%     y_log = log(met_rate);                     % log transform on y (m x 1)
%     beta_prime = X \ y_log;
%     R_sq = R_squared(X, y_log, beta_prime)
%     % Expected: R_sq close to 1 (tests path (h), log-log linearization)

    % Setup
    % Get the number of data points
    m = length(y);

    % Ensure y is a column vector (m x 1)
    y = y(:);

    % Create the ones vector (m x 1 column vector)
    j_vec = ones(m, 1);

    % Step 1: Compute SSE (sum of squared errors)
    % SSE = ||y - X * beta_prime||^2
    %   X * beta_prime: (m x n)(n x 1) = (m x 1)
    %   y - X * beta_prime: (m x 1) - (m x 1) = (m x 1)
    %   norm(...)^2: scalar
    SSE = norm(y - X * beta_prime)^2;

    % Step 2: Compute SST
    % Check if 1 is in CS(X) by testing if X can reproduce j_vec
    %   X \ j_vec:       (n x 1)  least-squares solution
    %   X * (X \ j_vec): (m x n)(n x 1) = (m x 1)  projection of j_vec onto CS(X)
    %   If the residual norm is near zero, then 1 is in CS(X)
    if (norm(X * (X \ j_vec) - j_vec) < 1e-12)
        % (h) 1 is in CS(X): use the centering matrix C = I - (1/m)*1*1'
        %   C:           (m x m) - (1/m)(m x 1)(1 x m) = (m x m)
        %   y' * C * y:  (1 x m)(m x m)(m x 1) = scalar
        C = eye(m) - (1/m) * (j_vec * j_vec');
        SST = y' * C * y;
    else
        % (i) 1 is NOT in CS(X): no centering needed
        %   y' * y: (1 x m)(m x 1) = scalar
        SST = y' * y;
    end

    % Step 3: Compute R^2
    % (j) Coefficient of determination
    %   R^2 = 1 - SSE / SST: scalar
    R_sq = 1 - SSE / SST;

end
