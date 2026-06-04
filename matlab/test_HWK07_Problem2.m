%% test_HWK07_Problem2.m
%  ================================================================
%  TEST SCRIPT FOR HWK_07_Problem_2_name.m
%  ================================================================
%  
%  HOW TO USE:
%    1. Place this file and HWK_07_Problem_2_name.m in the SAME folder
%    2. In MATLAB, cd to that folder
%    3. Run this script: >> test_HWK07_Problem2
%
%  This script runs 3 tests:
%    Test 1: Linear data        (degree 1) - should get R^2 ~ 1
%    Test 2: Quadratic data     (degree 2) - should get R^2 ~ 1
%    Test 3: Noisy cubic data   (degree 3) - should get R^2 close to 1
%  ================================================================

clear; clc; close all;

fprintf('============================================================\n');
fprintf('  TESTING HWK_07_Problem_2_name.m\n');
fprintf('============================================================\n\n');

%% ---- TEST 1: Perfect Linear Data ----
fprintf('>>> TEST 1: Linear data (degree 1)\n');
fprintf('    True model: y = 3 + 2x\n');

x1 = [1; 2; 3; 4; 5; 6; 7; 8];
y1 = 3 + 2*x1;   % perfect linear data, no noise

[beta1] = HWK_07_Problem_2_name(1, x1, y1);

fprintf('    Expected beta: [3.000000, 2.000000]\n');
fprintf('    Got beta:      [%.6f, %.6f]\n', beta1(1), beta1(2));

if abs(beta1(1) - 3) < 1e-6 && abs(beta1(2) - 2) < 1e-6
    fprintf('    RESULT: PASS\n\n');
else
    fprintf('    RESULT: FAIL\n\n');
end

%% ---- TEST 2: Quadratic Data with Small Noise ----
fprintf('>>> TEST 2: Quadratic data (degree 2)\n');
fprintf('    True model: y = 1 + 0.5x + 0.3x^2 + small noise\n');

x2 = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
y2 = 1 + 0.5*x2 + 0.3*x2.^2 + 0.1*randn(size(x2));

[beta2] = HWK_07_Problem_2_name(2, x2, y2);

fprintf('    Expected beta ~ [1.0, 0.5, 0.3]\n');
fprintf('    Got beta:       [%.6f, %.6f, %.6f]\n', beta2(1), beta2(2), beta2(3));
fprintf('    (Close but not exact due to noise - this is expected)\n');
fprintf('    RESULT: Check that R^2 is close to 1 in the figure title\n\n');

%% ---- TEST 3: Cubic Data ----
fprintf('>>> TEST 3: Cubic data (degree 3)\n');
fprintf('    True model: y = 2 - x + 0.5x^2 + 0.1x^3\n');

x3 = [0; 1; 2; 3; 4; 5; 6; 7; 8];
y3 = 2 - x3 + 0.5*x3.^2 + 0.1*x3.^3;   % perfect cubic

[beta3] = HWK_07_Problem_2_name(3, x3, y3);

fprintf('    Expected beta: [2.000000, -1.000000, 0.500000, 0.100000]\n');
fprintf('    Got beta:      [%.6f, %.6f, %.6f, %.6f]\n', ...
         beta3(1), beta3(2), beta3(3), beta3(4));

tol = 1e-6;
if abs(beta3(1)-2)<tol && abs(beta3(2)-(-1))<tol && ...
   abs(beta3(3)-0.5)<tol && abs(beta3(4)-0.1)<tol
    fprintf('    RESULT: PASS\n\n');
else
    fprintf('    RESULT: FAIL\n\n');
end

%% ---- SUMMARY ----
fprintf('============================================================\n');
fprintf('  ALL TESTS COMPLETE\n');
fprintf('  You should see 3 figures with data points and fitted curves.\n');
fprintf('  Each figure title shows the R^2 value.\n');
fprintf('============================================================\n');
fprintf('\n');
fprintf('  NEXT STEPS:\n');
fprintf('  1. Rename HWK_07_Problem_2_name.m to use your classmate''s name\n');
fprintf('     Example: HWK_07_Problem_2_Cortes.m\n');
fprintf('  2. Also rename the function line INSIDE the file to match\n');
fprintf('  3. Load your classmate''s .mat file:\n');
fprintf('     >> load(''classmate_file.mat'')\n');
fprintf('  4. Call the function with their data:\n');
fprintf('     >> [beta] = HWK_07_Problem_2_Cortes(poly_deg, x_vals_data, y_vals_data)\n');
fprintf('\n');
