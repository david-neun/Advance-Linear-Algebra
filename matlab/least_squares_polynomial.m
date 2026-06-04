% least_squares_polynomial.m
%
% Name:   David Cortes
% Date:   March 10, 2026
% Course: MATH 515 – Linear Algebra
%
% Description:
%   This script finds a least-squares polynomial fit of any specified
%   degree to a set of data points. It constructs a coefficient matrix
%   of powers of x and solves A*beta = y in the least-squares sense.
%
%   The model is:
%       y = beta_0 + beta_1*x + beta_2*x^2 + ... + beta_n*x^n
%
% Variables:
%   x_vals_data       – column vector (m x 1) of x-values of the data
%   y_vals_data       – column vector (m x 1) of y-values of the data
%   poly_deg          – scalar integer, degree of the fitting polynomial
%   num_points        – number of data points
%   coeff_matrix      – (m x (poly_deg+1)) matrix with columns x.^j
%   beta_vals         – column vector of least-squares coefficients

%% ---- Test Case 1: Linear fit (degree 1) --------------------------------
% Data that lie exactly on the line y = 2 + 3*x
% The least-squares solution should recover beta_0 = 2, beta_1 = 3.
fprintf('=== Test Case 1: Exact linear data (degree 1) ===\n');

x_test1 = [0; 1; 2; 3; 4];
y_test1 = [2; 5; 8; 11; 14];          % y = 2 + 3*x

beta_test1 = least_squares_fit(x_test1, y_test1, 1);
fprintf('Expected: beta_0 = 2, beta_1 = 3\n');
fprintf('Computed: beta_0 = %.6f, beta_1 = %.6f\n\n', ...
    beta_test1(1), beta_test1(2));

%% ---- Test Case 2: Quadratic fit (degree 2) -----------------------------
% Data that lie exactly on y = 1 - 2*x + 0.5*x^2
% The least-squares solution should recover those coefficients.
fprintf('=== Test Case 2: Exact quadratic data (degree 2) ===\n');

x_test2 = [-2; -1; 0; 1; 2; 3];
y_test2 = 1 - 2*x_test2 + 0.5*x_test2.^2;

beta_test2 = least_squares_fit(x_test2, y_test2, 2);
fprintf('Expected: beta_0 = 1, beta_1 = -2, beta_2 = 0.5\n');
fprintf('Computed: beta_0 = %.6f, beta_1 = %.6f, beta_2 = %.6f\n\n', ...
    beta_test2(1), beta_test2(2), beta_test2(3));

%% ---- Test Case 3: Constant fit (degree 0) -------------------------------
% Fitting a degree-0 polynomial is just computing the mean of the y-data.
fprintf('=== Test Case 3: Constant fit (degree 0) ===\n');

x_test3 = [1; 2; 3; 4; 5];
y_test3 = [4; 6; 8; 2; 10];           % mean = 6

beta_test3 = least_squares_fit(x_test3, y_test3, 0);
fprintf('Expected: beta_0 = %.6f (mean of y-data)\n', mean(y_test3));
fprintf('Computed: beta_0 = %.6f\n\n', beta_test3(1));

%% ---- Test Case 4: Cubic fit with noisy data (degree 3) -----------------
% Use data from y = 2 - x + 3*x^2 + 0.5*x^3, then add small noise.
% The fitted coefficients should be close to [2; -1; 3; 0.5].
fprintf('=== Test Case 4: Approximate cubic fit with noisy data (degree 3) ===\n');

x_test4 = (-3:0.5:3)';
y_exact4 = 2 - x_test4 + 3*x_test4.^2 + 0.5*x_test4.^3;
rng(42);                               % fix seed for reproducibility
y_test4 = y_exact4 + 0.1*randn(size(y_exact4));  % add small noise

beta_test4 = least_squares_fit(x_test4, y_test4, 3);
fprintf('Expected (approx): beta_0 ≈ 2, beta_1 ≈ -1, beta_2 ≈ 3, beta_3 ≈ 0.5\n');
fprintf('Computed: beta_0 = %.6f, beta_1 = %.6f, beta_2 = %.6f, beta_3 = %.6f\n\n', ...
    beta_test4(1), beta_test4(2), beta_test4(3), beta_test4(4));

%% ---- Test Case 5: Degree-11 polynomial fit ------------------------------
% Exact data from an 11th-degree polynomial with known coefficients.
% With enough data points, the least-squares solution should recover them.
fprintf('=== Test Case 5: Exact degree-11 polynomial ===\n');

% Use 20 evenly spaced points in [-1, 1] (more points than coefficients)
x_test5 = linspace(-1, 1, 20)';

% Known coefficients: beta_0 through beta_11
true_betas5 = [1; -2; 3; 0.5; -1; 0.25; 2; -0.3; 0.1; -0.05; 0.01; 0.004];

% Evaluate the polynomial at those x-values
y_test5 = zeros(size(x_test5));
for j = 0:11
    y_test5 = y_test5 + true_betas5(j+1) * x_test5.^j;
end

beta_test5 = least_squares_fit(x_test5, y_test5, 11);

fprintf('Expected vs Computed coefficients:\n');
for j = 0:11
    fprintf('  beta_%2d: expected = %10.6f, computed = %10.6f\n', ...
        j, true_betas5(j+1), beta_test5(j+1));
end
fprintf('\n');

%% ---- Application: Sixth-degree polynomial fit ---------------------------
fprintf('=== Sixth-degree polynomial fit to 415_sixth_degree_data ===\n');

% Load the data file (creates x_vals_data and y_vals_data)
load('415_sixth_degree_data');

% Fit a sixth-degree polynomial to the loaded data
beta_vals = least_squares_fit(x_vals_data, y_vals_data, 6);

% Display the resulting coefficients
fprintf('Least-squares coefficients (beta_0 through beta_6):\n');
for k = 0:6
    fprintf('  beta_%d = %14.8f\n', k, beta_vals(k+1));
end

%% ========================================================================
%  Local function: least_squares_fit
%  ========================================================================
function beta_vals = least_squares_fit(x_vals_data, y_vals_data, poly_deg)
% least_squares_fit  Compute the least-squares polynomial fit of any
%   specified degree to a set of data points.
%
% Inputs:
%   x_vals_data  – column vector (m x 1) of x-values
%   y_vals_data  – column vector (m x 1) of y-values
%   poly_deg     – degree of the polynomial
%
% Output:
%   beta_vals    – column vector ((poly_deg+1) x 1) of coefficients

% Ensure the input data are column vectors
x_vals_data = x_vals_data(:);
y_vals_data = y_vals_data(:);

% Number of data points
num_points = length(x_vals_data);

% Build the coefficient matrix A of size (num_points x (poly_deg + 1))
% Column j corresponds to x.^(j-1), for j = 1, 2, ..., poly_deg+1
coeff_matrix = zeros(num_points, poly_deg + 1);
for col = 0:poly_deg
    coeff_matrix(:, col + 1) = x_vals_data .^ col;
end

% Solve the least-squares problem A * beta = y using the backslash
% operator, which uses QR factorization internally and is more
% numerically stable than forming the normal equations (A'*A)\(A'*y)
beta_vals = coeff_matrix \ y_vals_data;
end
