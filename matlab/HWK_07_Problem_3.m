% HWK_07_Problem_3.m
%
% David Cortes
% March 30, 2026
%
% Description:
%   Problem 3: Cricket chirp rate vs temperature analysis.
%
%   Part (a): Fits polynomials of degrees 1 and 2 to the cricket data
%   with chirp rate as the independent variable (x) and temperature as
%   the dependent variable (y). Plots the data and fitted curves, reports
%   beta-values, R^2, and determines the appropriate polynomial degree.
%
%   Part (b): Compares the least-squares model to the rule of thumb:
%   "Add 40 to the number of chirps in 13 seconds to get temperature (F)."
%
% Variables:
%   cricket_data   - loaded data matrix (m x 2)
%   chirp_rate     - chirping rate in chirps/minute, independent variable (m x 1)
%   temperature    - temperature in degrees Fahrenheit, dependent variable (m x 1)
%   m              - number of data points (scalar)
%   poly_deg       - degree of the polynomial (scalar)
%   X              - Vandermonde design matrix (m x n)
%   beta_vals_1    - degree 1 least-squares coefficients (2 x 1)
%   beta_vals_2    - degree 2 least-squares coefficients (3 x 1)
%   R_sq_1         - R^2 for degree 1 fit (scalar)
%   R_sq_2         - R^2 for degree 2 fit (scalar)
%
% Dependencies:
%   least_squares_fit.m, R_squared.m, cricket_data.mat

clear; clc; close all;

%  Load data
load('cricket_data.mat');

% chirp rate is independent (x), temperature is dependent (y)
chirp_rate   = cricket_data(:, 2);   % chirps/minute (m x 1)
temperature  = cricket_data(:, 1);   % degrees Fahrenheit (m x 1)
m = length(chirp_rate);              % number of data points (scalar)

fprintf('  Problem 3: Cricket Chirp Rate vs Temperature\n');
fprintf('  Number of data points: %d\n', m);

%  Part (a): Fit polynomials of degree 1 and degree 2
fprintf('Part (a): Polynomial Fit\n\n');

% Degree 1 fit 
poly_deg_1 = 1;
n_1 = poly_deg_1 + 1;                              % n = 2 coefficients
beta_vals_1 = least_squares_fit(chirp_rate, temperature, poly_deg_1); % (2 x 1)

% Build Vandermonde matrix for R^2 (m x 2)
X_1 = zeros(m, n_1);
for j = 1:n_1
    X_1(:, j) = chirp_rate.^(j-1);
end
R_sq_1 = R_squared(X_1, temperature, beta_vals_1);  % scalar

fprintf('  Degree 1 polynomial:\n');
fprintf('    beta_0 = %.6f\n', beta_vals_1(1));
fprintf('    beta_1 = %.6f\n', beta_vals_1(2));
fprintf('    T(x) = %.6f + %.6f * x\n', beta_vals_1(1), beta_vals_1(2));
fprintf('    R^2 = %.6f\n\n', R_sq_1);

% Degree 2 fit 
poly_deg_2 = 2;
n_2 = poly_deg_2 + 1;                              % n = 3 coefficients
beta_vals_2 = least_squares_fit(chirp_rate, temperature, poly_deg_2); % (3 x 1)

% Build Vandermonde matrix for R^2 (m x 3)
X_2 = zeros(m, n_2);
for j = 1:n_2
    X_2(:, j) = chirp_rate.^(j-1);
end
R_sq_2 = R_squared(X_2, temperature, beta_vals_2);  % scalar

fprintf('  Degree 2 polynomial:\n');
fprintf('    beta_0 = %.6f\n', beta_vals_2(1));
fprintf('    beta_1 = %.6f\n', beta_vals_2(2));
fprintf('    beta_2 = %.6f\n', beta_vals_2(3));
fprintf('    T(x) = %.6f + %.6f * x + %.6f * x^2\n', ...
    beta_vals_2(1), beta_vals_2(2), beta_vals_2(3));
fprintf('    R^2 = %.6f\n\n', R_sq_2);

% Degree selection answer
fprintf('  Degree Selection:\n');
fprintf('    R^2 (degree 1) = %.6f\n', R_sq_1);
fprintf('    R^2 (degree 2) = %.6f\n', R_sq_2);
fprintf('    Improvement from degree 1 to 2: %.6f\n', R_sq_2 - R_sq_1);
fprintf('    A degree 1 polynomial captures the trend of the data.\n');
fprintf('    Degree 2 provides negligible improvement in R^2 and\n');
fprintf('    would overfit the data.\n\n');

% Smooth curves for plotting
x_plot = linspace(min(chirp_rate), max(chirp_rate), 200)'; % (200 x 1)

% Degree 1 smooth curve
X_plot_1 = [ones(200,1), x_plot];                          % (200 x 2)
y_plot_1 = X_plot_1 * beta_vals_1;                         % (200 x 1)

% Degree 2 smooth curve
X_plot_2 = [ones(200,1), x_plot, x_plot.^2];               % (200 x 3)
y_plot_2 = X_plot_2 * beta_vals_2;                         % (200 x 1)

% Plot: both fits on the same figure
figure;
plot(chirp_rate, temperature, 'ko', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'b', 'DisplayName', 'Data Points');
hold on;
plot(x_plot, y_plot_1, 'r-', 'LineWidth', 2, ...
     'DisplayName', sprintf('Degree 1 (R^2 = %.4f)', R_sq_1));
plot(x_plot, y_plot_2, 'g--', 'LineWidth', 2, ...
     'DisplayName', sprintf('Degree 2 (R^2 = %.4f)', R_sq_2));
hold off;
xlabel('Chirp Rate (chirps/minute)');
ylabel('Temperature (°F)');
title('Problem 3(a): Cricket Chirp Rate vs Temperature');
legend('Location', 'best');
grid on;


% Part (b): Rule of Thumb comparison
fprintf('Part (b): Rule of Thumb Comparison\n\n');

% Rule of thumb: temperature = (chirps in 13 seconds) + 40
% Since chirp_rate is in chirps/minute:
%   chirps in 13 seconds = chirp_rate * (13/60)
%   T_rule = 40 + (13/60) * chirp_rate
rule_slope     = 13/60;    % scalar
rule_intercept = 40;       % scalar

fprintf('  Rule of thumb: T = (chirps in 13 sec) + 40\n');
fprintf('  Converted to chirps/min: T = %.6f + %.6f * x\n\n', ...
    rule_intercept, rule_slope);

fprintf('  Model (degree 1):  T = %.6f + %.6f * x\n', ...
    beta_vals_1(1), beta_vals_1(2));
fprintf('  Rule of thumb:     T = %.6f + %.6f * x\n\n', ...
    rule_intercept, rule_slope);

fprintf('  Comparison:\n');
fprintf('    Intercept: Model = %.6f, Rule = %.6f, Diff = %.6f\n', ...
    beta_vals_1(1), rule_intercept, abs(beta_vals_1(1) - rule_intercept));
fprintf('    Slope:     Model = %.6f, Rule = %.6f, Diff = %.6f\n\n', ...
    beta_vals_1(2), rule_slope, abs(beta_vals_1(2) - rule_slope));

% --- Plot: model vs rule of thumb ---
y_rule_plot = rule_intercept + rule_slope * x_plot;  % (200 x 1)

figure;
plot(chirp_rate, temperature, 'ko', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'b', 'DisplayName', 'Data Points');
hold on;
plot(x_plot, y_plot_1, 'r-', 'LineWidth', 2, ...
     'DisplayName', sprintf('Model: T = %.2f + %.4f x', ...
     beta_vals_1(1), beta_vals_1(2)));
plot(x_plot, y_rule_plot, 'm--', 'LineWidth', 2, ...
     'DisplayName', sprintf('Rule of Thumb: T = 40 + %.4f x', rule_slope));
hold off;
xlabel('Chirp Rate (chirps/minute)');
ylabel('Temperature (°F)');
title('Problem 3(b): Model vs Rule of Thumb');
legend('Location', 'best');
grid on;

fprintf('  Problem 3 Complete. Check the two figures.\n');