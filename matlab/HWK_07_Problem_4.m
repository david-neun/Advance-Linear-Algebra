% HWK_07_Problem_4.m
%
% David Cortes
% March 30, 2026
%
% Description:
%   Problem 4: Metabolic rate vs mass analysis.
%
%   The data does no fit a polynomial, but can be made linear using
%   logarithmic transforms. This script tries three transforms:
%     Transform 1: y = beta_0 + beta_1 * log(x)        (log on x only)
%     Transform 2: log(y) = beta_0 + beta_1 * x        (log on y only)
%     Transform 3: log(y) = beta_0 + beta_1 * log(x)   (log on both)
%
%   Fits a line to each, compares R^2 to find the best transform, then
%   derives and plots y = f(x) for the non-transformed data.
%
% Variables:
%   metabolic_rate_data - loaded data matrix (m x 2)
%   mass                - mass in kg, independent variable (m x 1)
%   met_rate            - metabolic rate in cal/day, dependent variable (m x 1)
%   m                   - number of data points (scalar)
%   X_t1, X_t2, X_t3   - design matrics for each transform (m x 2)
%   beta_t1, beta_t2, beta_t3 - coefficients for each transform (2 x 1)
%   R_sq_t1, R_sq_t2, R_sq_t3 - R^2 for each transform (scalar)
%
% Dependencies:
%   R_squared.m, metabolic_rate_data.mat

clear; clc; close all;

%  Load data
load('metabolic_rate_data.mat');

mass     = metabolic_rate_data(:, 1);   % mass in kg (m x 1)
met_rate = metabolic_rate_data(:, 2);   % metabolic rate in cal/day (m x 1)
m = length(mass);                       % number of data points (scalar)

fprintf('  Problem 4: Metabolic Rate vs Mass\n');
fprintf('  Number of data points: %d\n', m);

%  Transform 1: y = beta_0 + beta_1 * log(x)
%  Only the independent variable is transformed
%  Design matrix: X = [1, log(x)]
X_t1 = [ones(m, 1), log(mass)];         % (m x 2)
beta_t1 = X_t1 \ met_rate;              % (2 x 1)
R_sq_t1 = R_squared(X_t1, met_rate, beta_t1);  % scalar

fprintf('Transform 1: y = beta_0 + beta_1 * log(x)\n');
fprintf('  beta_0 = %.6f\n', beta_t1(1));
fprintf('  beta_1 = %.6f\n', beta_t1(2));
fprintf('  R^2 = %.6f\n\n', R_sq_t1);

%  Transform 2: log(y) = beta_0 + beta_1 * x
%  Only the dependent variable is transformed
%  design matrix: X = [1, x],  y_transformed = log(y)
X_t2 = [ones(m, 1), mass];              % (m x 2)
log_met_rate = log(met_rate);            % (m x 1)
beta_t2 = X_t2 \ log_met_rate;          % (2 x 1)
R_sq_t2 = R_squared(X_t2, log_met_rate, beta_t2);  % scalar

fprintf('Transform 2: log(y) = beta_0 + beta_1 * x\n');
fprintf('  beta_0 = %.6f\n', beta_t2(1));
fprintf('  beta_1 = %.6f\n', beta_t2(2));
fprintf('  R^2 = %.6f\n\n', R_sq_t2);

%  Transform 3: log(y) = beta_0 + beta_1 * log(x)
%  both variables are transformed
%  design matrix: X = [1, log(x)],  y_transformed = log(y)
X_t3 = [ones(m, 1), log(mass)];         % (m x 2)
beta_t3 = X_t3 \ log_met_rate;          % (2 x 1)
R_sq_t3 = R_squared(X_t3, log_met_rate, beta_t3);  % scalar

fprintf('Transform 3: log(y) = beta_0 + beta_1 * log(x)\n');
fprintf('  beta_0 = %.6f\n', beta_t3(1));
fprintf('  beta_1 = %.6f\n', beta_t3(2));
fprintf('  R^2 = %.6f\n\n', R_sq_t3);

%  Determine best transform
fprintf(' Transform Comparison \n');
fprintf('  Transform 1 (log x only):   R^2 = %.6f\n', R_sq_t1);
fprintf('  Transform 2 (log y only):   R^2 = %.6f\n', R_sq_t2);
fprintf('  Transform 3 (log x, log y): R^2 = %.6f\n\n', R_sq_t3);

% ========================================================================
%  Plot the three transforms to verify linearity
% ========================================================================

% Figure 1: Transform 1 - y vs log(x)
figure;
plot(log(mass), met_rate, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
x_line = linspace(min(log(mass)), max(log(mass)), 200)';  % (200 x 1)
y_line = beta_t1(1) + beta_t1(2) * x_line;                % (200 x 1)
plot(x_line, y_line, 'r-', 'LineWidth', 2);
hold off;
xlabel('log(mass)');
ylabel('Metabolic Rate (cal/day)');
title(sprintf('Transform 1: y vs log(x), R^2 = %.4f', R_sq_t1));
grid on;

% Figure 2: Transform 2 - log(y) vs x
figure;
plot(mass, log_met_rate, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
x_line2 = linspace(min(mass), max(mass), 200)';            % (200 x 1)
y_line2 = beta_t2(1) + beta_t2(2) * x_line2;               % (200 x 1)
plot(x_line2, y_line2, 'r-', 'LineWidth', 2);
hold off;
xlabel('Mass (kg)');
ylabel('log(Metabolic Rate)');
title(sprintf('Transform 2: log(y) vs x, R^2 = %.4f', R_sq_t2));
grid on;

% Figure 3: Transform 3 - log(y) vs log(x)
figure;
plot(log(mass), log_met_rate, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
y_line3 = beta_t3(1) + beta_t3(2) * x_line;                % (200 x 1)
plot(x_line, y_line3, 'r-', 'LineWidth', 2);
hold off;
xlabel('log(mass)');
ylabel('log(Metabolic Rate)');
title(sprintf('Transform 3: log(y) vs log(x), R^2 = %.4f', R_sq_t3));
grid on;

% Non-transformed function y = f(x)
% Transform 3 gives: log(y) = beta_0 + beta_1 * log(x)
%   => y = exp(beta_0 + beta_1 * log(x))
%   => y = exp(beta_0) * x^beta_1
%
% This is a power law: y = a * x^b
%   where a = exp(beta_0) and b = beta_1

a = exp(beta_t3(1));   % scalar
b = beta_t3(2);        % scalar

fprintf(' Non-Transformed Function \n');
fprintf('  From Transform 3: log(y) = %.6f + %.6f * log(x)\n', ...
    beta_t3(1), beta_t3(2));
fprintf('  Exponentiate both sides:\n');
fprintf('    y = exp(%.6f) * x^(%.6f)\n', beta_t3(1), beta_t3(2));
fprintf('    y = %.6f * x^(%.6f)\n\n', a, b);

% Figure 4: Non-transformed data with f(x) overlay
figure;
plot(mass, met_rate, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'b', ...
     'DisplayName', 'Data Points');
hold on;
x_smooth = linspace(min(mass), max(mass), 200)';  % (200 x 1)
y_fit = a * x_smooth.^b;                          % (200 x 1)
plot(x_smooth, y_fit, 'r-', 'LineWidth', 2, ...
     'DisplayName', sprintf('y = %.4f * x^{%.4f}', a, b));
hold off;
xlabel('Mass (kg)');
ylabel('Metabolic Rate (cal/day)');
title('Problem 4: Metabolic Rate vs Mass (Non-Transformed)');
legend('Location', 'best');
grid on;

fprintf(' Problem 4 Complete. Check the four figures.\n');
