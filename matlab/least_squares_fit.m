function beta_vals = least_squares_fit(x_vals_data, y_vals_data, poly_deg)
% David Cortes
% March 10, 2026
%
% least_squares_fit - Compute the least-squares polynomial fit of any
% specified degree to a set of data points.
%
% Constructs a Vandermonde matrix and solves A*beta = y in the
% least-squares sense using MATLAB's backslash operator.
%
% USAGE:
%   beta_vals = least_squares_fit(x_vals_data, y_vals_data, poly_deg)
%
% INPUTS:
%   x_vals_data  - column vector (m x 1) of x-values
%   y_vals_data  - column vector (m x 1) of y-values
%   poly_deg     - degree of the polynomial (scalar)
%
% OUTPUT:
%   beta_vals    - column vector ((poly_deg+1) x 1) of coefficients
%
% TEST CASES:
%
%   Test 1: Exact linear data (degree 1)
%     x = [0; 1; 2; 3; 4];
%     y = [2; 5; 8; 11; 14];
%     beta = least_squares_fit(x, y, 1)
%     % Expected: beta_0 = 2, beta_1 = 3
%
%   Test 2: Exact quadratic data (degree 2)
%     x = [-2; -1; 0; 1; 2; 3];
%     y = 1 - 2*x + 0.5*x.^2;
%     beta = least_squares_fit(x, y, 2)
%     % Expected: beta_0 = 1, beta_1 = -2, beta_2 = 0.5

    % Ensure the input data are column vectors (m x 1)
    x_vals_data = x_vals_data(:);
    y_vals_data = y_vals_data(:);

    % Number of data points
    num_points = length(x_vals_data);

    % Build the Vandermonde matrix A of size (num_points x (poly_deg + 1))
    % Column j corresponds to x.^(j-1), for j = 1, 2, ..., poly_deg+1
    coeff_matrix = zeros(num_points, poly_deg + 1);
    for col = 0:poly_deg
        coeff_matrix(:, col + 1) = x_vals_data .^ col;
    end

    % Solve the least-squares problem A * beta = y using backslash
    beta_vals = coeff_matrix \ y_vals_data;

end
