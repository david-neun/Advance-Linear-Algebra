% test_Gram_Schmidt.m

A = [ 1  3  5;
     -1 -3  1;
      0  2  3;
      1  5  2;
      1  5  8 ];

[Q, R] = Gram_Schmidt(A);

% --- Basic checks ---
recon_err = norm(A - Q*R, 'fro');      % should be ~ 0
orth_err  = norm(Q'*Q - eye(size(Q,2)), 'fro');  % should be ~ 0

% Check R is (numerically) upper triangular
lower_part = tril(R, -1);              % entries below diagonal
tri_err = norm(lower_part, 'fro');     % should be ~ 0

fprintf('||A - Q*R||_F     = %.3e\n', recon_err);
fprintf('||Q''*Q - I||_F    = %.3e\n', orth_err);
fprintf('||tril(R,-1)||_F   = %.3e\n', tri_err);

% --- Optional: compare with MATLAB QR (economy) ---
[Qm, Rm] = qr(A, 0);

fprintf('\nCompare with MATLAB qr(A,0):\n');
fprintf('||A - Qm*Rm||_F    = %.3e\n', norm(A - Qm*Rm, 'fro'));
fprintf('||Qm''*Qm - I||_F   = %.3e\n', norm(Qm'*Qm - eye(size(Qm,2)), 'fro'));

% Note: Q and R are not unique (sign differences possible).
% A quick way to compare is the reconstruction error above.
