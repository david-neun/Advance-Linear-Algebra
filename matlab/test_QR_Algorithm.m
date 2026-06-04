% test_QR_Algorithm.m

A0 = [ 4  1  2;
       1  3  0;
       2  0  2 ];      % symmetric 3x3

numIter = 40;

[Aend, eigApprox, A_hist] = QR_Algorithm(A0, numIter);

fprintf('Final A_%d:\n', numIter);
disp(Aend);

fprintf('Diagonal (eigenvalue approx):\n');
disp(eigApprox);

fprintf('Compare with MATLAB eig(A0):\n');
disp(eig(A0));

% Optional: check it's "almost" upper triangular at the end
fprintf('||tril(Aend,-1)||_F = %.3e\n', norm(tril(Aend,-1), 'fro'));
