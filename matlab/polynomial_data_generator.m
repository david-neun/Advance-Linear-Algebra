function [poly_deg, x_vals_data, y_vals_data, c_vals] = polynomial_data_generator()
% polynomial_data_generator  Generate synthetic polynomial data with noise.
%
% Name:   David Cortes
% Date:   March 12, 2026
% Course: MATH 515 – Linear Algebra
%
% Description:
%   This function generates a data set that approximately follows a
%   polynomial of a chosen degree. It picks d+1 exact points, solves
%   for the polynomial coefficients using a linear system, then
%   generates n random points on that polynomial and adds normally
%   distributed noise. Finally it plots the result.
%
% Outputs:
%   poly_deg    – degree of the polynomial (integer, 4 <= d <= 10)
%   x_vals_data – column vector (n x 1) of randomly generated x-values
%   y_vals_data – column vector (n x 1) of y-values with added noise
%   c_vals      – column vector ((poly_deg+1) x 1) of polynomial
%                 coefficients [c_0; c_1; ... ; c_d]
%
% Usage:
%   [poly_deg, x_vals_data, y_vals_data, c_vals] = polynomial_data_generator()

%% (a) Choose a degree for the polynomial (4 <= d <= 10)
poly_deg = 5;

%% (b) Choose d+1 points and compute coefficients via a linear system
% Pick d+1 equally spaced x-values on [0, 10] as the exact data points
num_exact = poly_deg + 1;
x_exact = linspace(0, 10, num_exact)';

% Choose y-values for those points (arbitrary choices)
y_exact = [3; -1; 2; 5; 0; 4];

% Build the coefficient matrix for the linear system:
%   c_0 + c_1*x + c_2*x^2 + ... + c_d*x^d = y
% Each row corresponds to one exact data point
exact_matrix = zeros(num_exact, num_exact);
for col = 0:poly_deg
    exact_matrix(:, col + 1) = x_exact .^ col;
end

% Solve the linear system (exact fit, not least squares)
c_vals = exact_matrix \ y_exact;

%% (c) Choose the number of random data points (50 <= n <= 100)
num_points = 75;

%% (d) Randomly generate n x-values and evaluate the polynomial
% Seed the random number generator for reproducibility across runs
% using 'clock', I got this warning: 'clock' is not recommended. 
% With appropriate code changes, use 'datetime("now")! instead.
% I avoided to changed it.
rng(sum(clock), 'twister');

% Generate random x-values on the interval [0, 10]
a = 0;
b = 10;
x_vals_data = a + (b - a) * rand(num_points, 1);

% Evaluate the polynomial at those x-values using polyval
% polyval expects coefficients in descending order, so flip c_vals
y_clean = polyval(flip(c_vals), x_vals_data);

%% (e) Add normally distributed noise with mean 0 and small std dev
% Note: normrnd requires the Statistics and Machine Learning Toolbox.
% I installed it via MATLAB Add-Ons to run this code.
%
% Choose sigma proportional to the range of the clean y-values so the
% noise is large enough to be visible but small enough that the data
% still follows the shape of the polynomial.
y_range = max(y_clean) - min(y_clean);
sigma = 0.05 * y_range;
noise = normrnd(0, sigma, num_points, 1);
y_vals_data = y_clean + noise;

%% (f) Plot the noisy data points against the smooth polynomial curve
figure;

% Create a smooth curve for the polynomial
x_smooth = linspace(a, b, 500)';
y_smooth = polyval(flip(c_vals), x_smooth);

% Plot the smooth polynomial and the noisy data
plot(x_smooth, y_smooth, 'b-', 'LineWidth', 1.5); hold on;
plot(x_vals_data, y_vals_data, 'ro', 'MarkerSize', 5);
hold off;

xlabel('x');
ylabel('y');
title(sprintf('Degree-%d Polynomial with Noisy Data (\\sigma = %.2f)', ...
    poly_deg, sigma));
legend('True polynomial', 'Noisy data', 'Location', 'best');
grid on;

%% (g) Save x_vals_data, y_vals_data, and poly_deg to a .mat file
save('Cortes_data.mat', 'x_vals_data', 'y_vals_data', 'poly_deg');
fprintf('Data saved to Cortes_data.mat\n');

end
