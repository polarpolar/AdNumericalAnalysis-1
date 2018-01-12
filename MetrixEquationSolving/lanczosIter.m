function [ iter, x_res, error, intermediates ] = lanczosIter( M, b, err_thr, max_iter )
%LANCZOSITER this is the Lanczos solution of equation solving.
%   The Lanczos solution of metrix equation solving. The lanczos solution
%   is based on the Galerkin principle by searching equation solutions in
%   sub-spaces. In lanczos solution, we construct the sub-spaces with
%   krylov sub-spaces and the same search/verification spaces V = W. 
%   
%   In lanczos solution, we search sub-spaces from 1 to Msize dimensions.
%   Each space is constructed by lanczos process to solve space bases and
%   simplified equivalent equations. The termination condition is defined
%   by err_thr of Mx = b. To avoid in-convergency, a maximum iteration
%   number is defined. If not defined, the maximum iteration number is
%   limited to Msize, because the Msize-th metrix is retreat to the
%   original metrix, M.
%
%   After all, we give the iteration times of lanczos solution, the result,
%   the error and all of the intermediates variables.

% check arguments
[Msize, ~] = size(M);
if nargin <= 3; max_iter = Msize + 1; end

%% Lanczos Method to solve Equation
% accurate result
x_acc = M \ b;

% variables
x_iter = zeros(Msize, max_iter);
y_iter = zeros(Msize, max_iter);
r_iter = zeros(1, max_iter);
err2_iter = zeros(1, max_iter);
errA_iter = zeros(1, max_iter);

% initial value (offset for linear flow pattern)
%   Different from the iteration approachs using Ritz's Theory, Galerkin
%   Theory does not need to initialize a start point for iteration but
%   requires to modify the x_iter-s to a linear flow pattern by per-define
%   an offset (x_init = 0).
x_init = zeros(Msize, 1);
beta_init = norm(b - M * x_init);

%   Give the complete lanczos process. All the intermediate metrics are
%   defined as the k-th sub-metrix of Q.
[ Q, T ] = lanczosProc( M, b, Msize );
for iter = 1 : max_iter
    % Find the base metrix of the k-th sub-space. The base metrix is
    % constructed by the k-th sub-metrix of lanczos metrix Q, and T. 
    % get Q^(k) & T^(k) with Lanczos Process
    Q_iter = Q(:, 1:iter);
    T_iter = T(1:iter, 1:iter);
    if rank(T_iter) < iter
        % if the metrix is singular, we just need to move on to the next
        % metrix before reaching the largest.
        continue;
    end
    
    % Now we solve the simplified equation T * y = (mod(r0), 0, ..., 0)'
    e1 = zeros(iter, 1);    e1(1) = beta_init;
    y_iter(1:iter, iter) = T_iter \ e1;
    x_iter(:, iter) = x_init + Q_iter * y_iter(1:iter, iter);
    
    % Termination condiftion according to r_k
    %   In lanczos, r_k = -beta_k * q(k+1)' * (ek'*y_k);
    %   To simplified the computation:
    %                   |r_k| = |beta_k*(ek'*y_k)|
    if iter == max_iter
        beta_k = 0;
    else
        beta_k = T(iter, iter+1);
    end
    ek = zeros(iter, 1);
    ek(iter) = 1;
    r_iter(iter) = abs( beta_k*ek'*y_iter(1:iter, iter) );
    % Forward Error
    err2_iter(iter) = norm( x_iter(:, iter) - x_acc );
    errA_iter(iter) = Anorm( x_iter(:, iter) - x_acc, M );
    if r_iter(iter) < err_thr
        break
    end
end

x_res = x_iter(:, iter);
% this is error 2
error = r_iter(iter);

intermediates.iter_num = iter;
intermediates.x_iter = x_iter(:, 1:iter);
intermediates.r_iter = r_iter(1:iter)';
intermediates.err2_iter = err2_iter(1:iter)';
intermediates.errA_iter = errA_iter(1:iter)';

end

