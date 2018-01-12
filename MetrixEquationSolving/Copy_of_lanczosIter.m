function [ iter_num, x_res, error, intermediates ] = lanczosIter( M, b, err_thr, max_iter )
%LANCZOSITER this is the lanczos solution of equation solving.
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

% iterate variables
x_iter = zeros(Msize, max_iter);
y_iter = zeros(Msize, max_iter);
r_iter = zeros(1, max_iter);
err2_iter = zeros(1, max_iter);
errA_iter = zeros(1, max_iter);

% initial value
x_iter(:, 1) = zeros(Msize, 1);
r_iter(:, 1) = norm( b - M * x_iter(:, 1) );

%   Give the complete lanczos process. All the intermediate metrics are
%   defined as the k-th sub-metrix of Q.
[ Q, T ] = lanczosProc( M, b, Msize );
for iter_num = 2 : max_iter
    k = iter_num - 1;
    
    % Find the base metrix of the k-th sub-space. The base metrix is
    % constructed by the k-th sub-metrix of lanczos metrix Q, and T. 
    Q_iter = Q(:, 1:k);
    T_iter = T(1:k, 1:k);
    if rank(T_iter) < k
        % if the metrix is singular, we just need to move on to the next
        % metrix before reaching the largest.
        continue
    end
    
    % Now we solve the simplified equation T * y = (mod(r0), 0, ..., 0)'
    e1 = zeros(k, 1);
    e1(1) = r_iter(:, 1);
    y_iter(1:k, iter_num) = T_iter \ e1;
    x_iter(:, iter_num) = x_iter(:, 1) + Q_iter * y_iter(1:k, iter_num);
    
    % Calculate error
    %r = norm(M * Q_iter * y_iter(1:k, iter_num) - r_iter(:, 1));
    beta_k = T(k, k+1);
    ek = zeros(k, 1);
    ek(k) = 1;
    err2_iter(iter_num) = norm( x_iter(:, iter_num) - x_acc );
    errA_iter(iter_num) = Anorm( x_iter(:, iter_num) - x_acc, M );
    
    % termination condiftion according to r_k
    r_iter(iter_num) = abs( beta_k*ek'*y_iter(1:k, iter_num) );
    if r_iter(iter_num) < err_thr
        break
    end
end

error = err2_iter(iter_num);
x_res = x_iter(:, iter_num);
intermediates.iter_num = iter_num - 1;
intermediates.x_iter = x_iter(2:iter_num);
intermediates.r_iter = r_iter(1:iter_num);
intermediates.err2_iter = err2_iter(2:iter_num);
intermediates.errA_iter = errA_iter(2:iter_num);

% remove initial index
iter_num = iter_num - 1;
end

