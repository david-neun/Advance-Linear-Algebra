function [beta_vals] = HWK_07_Problem_2_name(poly_deg, x_vals_data, y_vals_data)
% HWK_07_Problem_2_name
%
% Fits a polynomial of degree poly_deg to the given data using
% least-squares, computes R^2, and plots the data with the fitted curve.
%
% USAGE:
%   [beta_vals] = HWK_07_Problem_2_name(poly_deg, x_vals_data, y_vals_data)
%
% INPUTS:
%   poly_deg    - degree of the polynomial to fit (scalar)
%   x_vals_data - independent variable data (m x 1)
%   y_vals_data - dependent variable data   (m x 1)
%
% OUTPUT:
%   beta_vals   - least-squares coefficients ((poly_deg+1) x 1)
%
% RENAME THIS FUNCTION AND FILE: Replace "name" with your classmate's
% name. For example: HWK_07_Problem_2_Cortes.m
% Remember: the function name MUST match the file name.

    % ================================================================
    % Step 0: Setup and dimension extraction
    % ================================================================
    m = length(y_vals_data);        % m = number of data points (scalar)
    n = poly_deg + 1;               % n = number of columns in X (scalar)

    % Ensure column vectors: (m x 1)
    x_vals_data = x_vals_data(:);   % force (m x 1)
    y_vals_data = y_vals_data(:);   % force (m x 1)

    % ================================================================
    % Step A: Build the Vandermonde design matrix X
    % ================================================================
    % X is (m x n), where n = poly_deg + 1
    % Column j of X is x_vals_data.^(j-1), for j = 1, ..., n
    %
    % Dimension check:
    %   Each column x_vals_data.^(j-1) is (m x 1)
    %   We stack n columns side by side -> X is (m x n)

    X = zeros(m, n);                         % pre-allocate (m x n)
    for j = 1:n
        X(:, j) = x_vals_data.^(j-1);       % (m x 1) into column j
    end
    % X is now (m x n) where n = poly_deg + 1

    % ================================================================
    % Step B: Solve the least-squares problem X * beta = y
    % ================================================================
    % MATLAB backslash solves the normal equations:
    %   (X^T X) beta = X^T y
    %   (n x m)(m x n) beta = (n x m)(m x 1)
    %   (n x n) beta = (n x 1)
    % Result: beta_vals is (n x 1)

    beta_vals = X \ y_vals_data;             % (n x 1)

    % ================================================================
    % Step C: Compute R^2 (Coefficient of Determination)
    % ================================================================

    % --- Compute SSE ---
    %   X * beta_vals: (m x n)(n x 1) = (m x 1)
    %   y - X*beta:    (m x 1) - (m x 1) = (m x 1)
    %   norm(...)^2:   scalar

    y_pred = X * beta_vals;                  % (m x 1)
    SSE = norm(y_vals_data - y_pred)^2;      % scalar

    % --- Determine SST ---
    % Check if 1 is in CS(X) [Problem 1 part (g)]
    %   j_vec:           (m x 1)
    %   X \ j_vec:       (n x 1)
    %   X*(X\j_vec):     (m x n)(n x 1) = (m x 1)
    %   residual:        (m x 1) - (m x 1) = (m x 1)

    j_vec = ones(m, 1);                      % (m x 1)

    if (norm(X * (X \ j_vec) - j_vec) < 1e-12)
        % 1 is in CS(X): use centering matrix
        % C = I - (1/m)*1*1^T:  (m x m)
        % SST = y^T * C * y:   (1xm)(mxm)(mx1) = scalar
        C = eye(m) - (1/m) * (j_vec * j_vec');
        SST = y_vals_data' * C * y_vals_data;
    else
        % 1 is NOT in CS(X): no centering
        % SST = y^T * y:  (1xm)(mx1) = scalar
        SST = y_vals_data' * y_vals_data;
    end

    % --- Compute R^2 ---
    R_squared = 1 - SSE / SST;               % scalar

    % ================================================================
    % Step D: Generate smooth curve for plotting
    % ================================================================
    x_plot = linspace(min(x_vals_data), max(x_vals_data), 200)';
    X_plot = zeros(length(x_plot), n);
    for j = 1:n
        X_plot(:, j) = x_plot.^(j-1);
    end
    y_plot = X_plot * beta_vals;             % (200 x 1)

    % ================================================================
    % Step E: Plot the data and the fitted polynomial
    % ================================================================
    figure;
    plot(x_vals_data, y_vals_data, 'ko', 'MarkerSize', 8, ...
         'MarkerFaceColor', 'b', 'DisplayName', 'Data Points');
    hold on;
    plot(x_plot, y_plot, 'r-', 'LineWidth', 2, ...
         'DisplayName', sprintf('Degree %d Polynomial', poly_deg));
    hold off;

    xlabel('x');
    ylabel('y');
    title(sprintf('Least-Squares Polynomial Fit (degree %d), R^2 = %.6f', ...
           poly_deg, R_squared));
    legend('Location', 'best');
    grid on;

    % ================================================================
    % Step F: Display results to the command window
    % ================================================================
    fprintf('\n--- Least-Squares Polynomial Fit Results ---\n');
    fprintf('Polynomial degree: %d\n', poly_deg);
    fprintf('Number of data points (m): %d\n', m);
    fprintf('Number of parameters  (n): %d\n', n);
    fprintf('\nDimension Summary:\n');
    fprintf('  X:          (%d x %d)\n', m, n);
    fprintf('  beta_vals:  (%d x 1)\n', n);
    fprintf('  y_pred:     (%d x 1)\n', m);
    fprintf('\nBeta values:\n');
    for j = 1:n
        fprintf('  beta_%d = %12.6f\n', j-1, beta_vals(j));
    end
    fprintf('\nPolynomial: y = ');
    for j = 1:n
        if j == 1
            fprintf('%.6f', beta_vals(j));
        else
            if beta_vals(j) >= 0
                fprintf(' + %.6f*x^%d', beta_vals(j), j-1);
            else
                fprintf(' - %.6f*x^%d', abs(beta_vals(j)), j-1);
            end
        end
    end
    fprintf('\n');
    fprintf('\nSSE = %.6f\n', SSE);
    fprintf('SST = %.6f\n', SST);
    fprintf('R^2 = %.6f\n', R_squared);
    fprintf('--------------------------------------------\n\n');

end
