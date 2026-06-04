function [beta_vals] = HWK_07_Problem_2_Hanson(poly_deg, x_vals_data, y_vals_data)
% David Cortes
% March 30, 2026
%
% HWK_07_Problem_2_Hanson
%
% Fits a polynomial of degree poly_deg to Hanson's data using
% least-squares, computes R^2, and plots the data with the fitted curve.
% Calls least_squares_fit.m and R_squared.m from previous work.
%
% USAGE:
%   [beta_vals] = HWK_07_Problem_2_Hanson(poly_deg, x_vals_data, y_vals_data)
%
% INPUTS:
%   poly_deg    - degree of the polynomial to fit (scalar)
%   x_vals_data - independent variable data (m x 1 column vector)
%   y_vals_data - dependent variable data   (m x 1 column vector)
%
% OUTPUT:
%   beta_vals   - least-squares coefficients ((poly_deg+1) x 1 column vector)
%
% TEST CASES:
%
%   Test 1: Using Hanson's data
%     load('Hanson_data.mat');
%     [beta_vals] = HWK_07_Problem_2_Hanson(poly_deg, x_vals_data, y_vals_data)
%
%   Test 2: Known quadratic data (degree 2, R^2 should equal 1)
%     x = [0; 1; 2; 3; 4];
%     y = [1; 2; 5; 10; 17];
%     [beta_vals] = HWK_07_Problem_2_Hanson(2, x, y)
%     % Expected: beta_0 = 1, beta_1 = 0, beta_2 = 1, R^2 = 1.0

    % ================================================================
    % Step 0: Setup
    % ================================================================
    m = length(y_vals_data);        % m = number of data points (scalar)
    n = poly_deg + 1;               % n = number of coefficients (scalar)

    % Ensure column vectors (m x 1)
    x_vals_data = x_vals_data(:);
    y_vals_data = y_vals_data(:);

    % ================================================================
    % Step A: Solve least-squares using least_squares_fit
    % ================================================================
    % Calls least_squares_fit.m which builds the Vandermonde matrix
    % and solves the system via backslash
    beta_vals = least_squares_fit(x_vals_data, y_vals_data, poly_deg);  % (n x 1)

    % ================================================================
    % Step B: Build X for R^2 computation
    % ================================================================
    % Reconstruct the Vandermonde matrix X (m x n) for R_squared
    X = zeros(m, n);
    for j = 1:n
        X(:, j) = x_vals_data.^(j-1);       % (m x 1) into column j
    end

    % ================================================================
    % Step C: Compute R^2 using R_squared function
    % ================================================================
    R_sq = R_squared(X, y_vals_data, beta_vals);  % scalar

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
    title(sprintf('Hanson Data: Degree %d Fit, R^2 = %.6f', ...
           poly_deg, R_sq));
    legend('Location', 'best');
    grid on;

    % ================================================================
    % Step F: Display results to the command window
    % ================================================================
    fprintf('\n--- Hanson Data: Least-Squares Polynomial Fit ---\n');
    fprintf('Polynomial degree: %d\n', poly_deg);
    fprintf('Number of data points (m): %d\n', m);
    fprintf('Number of parameters  (n): %d\n', n);
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
    fprintf('\nR^2 = %.6f\n', R_sq);
    fprintf('------------------------------------------------\n\n');

end
