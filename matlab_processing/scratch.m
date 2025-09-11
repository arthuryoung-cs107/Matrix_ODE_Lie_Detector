clear

A = rand(30,10);
b = rand(30,1);

[U,S,V] = svd(A,'econ'); % thin SVD of A

%{
    Since the entries of A are random, A is almost surely full rank.
    Therefore, using S as in the thin SVD of A, S^-1 is well defined.
%}
A_dagg = V *( inv(S) )*U'; % by definition of the Moore--Penrose pseudoinverse

x_hat = A_dagg*b;
rsqrd =  sum( (A*x_hat-b).^2 );

eps_scl = 1e-2; % scale for small perturbation of x_hat
d1 = eps_scl*rand(10,1);
d2 = eps_scl*rand(10,1);
d3 = eps_scl*rand(10,1);

rsqrd1 = sum( (A*(x_hat + d1)-b).^2 );
rsqrd2 = sum( (A*(x_hat + d2)-b).^2 );
rsqrd3 = sum( (A*(x_hat + d3)-b).^2 );

fprintf('Minimum residual is r = %.3e. r1/r = %.3f, r2/r = %.3f, r2/r = %.3f.\n', ...
rsqrd, ...
rsqrd1 / rsqrd, ...
rsqrd2 / rsqrd, ...
rsqrd3 / rsqrd );
