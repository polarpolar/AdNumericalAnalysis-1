function [ iter, x_res, error, intermediates ] = cgIter( M, b, err_thr, max_iter )
%CGITER Summary of this function goes here
%   Detailed explanation goes here

%% C-G
% accurate value
x_acc = M \ b;
[Msize, ~] = size(M);

% variables
x_iter = zeros(Msize, max_iter);
r_iter = zeros(Msize, max_iter);
rn_iter = zeros(1, max_iter);
p_iter = zeros(Msize, max_iter);
err2_iter = zeros(1, max_iter);
errA_iter = zeros(1, max_iter);

% initial value
x_iter(:, 1) = zeros(Msize, 1);
r_iter(:, 1) = b - M * x_iter(:, 1);
p_iter(:, 1) = r_iter(:, 1);

% iteration process
for iter = 1:max_iter
    % Termination Condition, the error is limited by thr.
%   err2_iter(iter_num) = norm(r_iter(:, iter_num));
% 	errA_iter(iter_num) = sqrt(Anorm(r_iter(:, iter_num), M));
    err2_iter(iter) = norm(x_iter(:, iter) - x_acc);
    errA_iter(iter) = sqrt(Anorm(x_iter(:, iter) - x_acc, M));
    rn_iter(iter) = norm(r_iter(:, iter));
    if rn_iter(iter) < err_thr
        break
    end
    
    % Construct x(:, iter+1) with base vector p(:, iter)
    alpha_iter = norm(r_iter(:, iter))^2 / Anorm(p_iter(:, iter), M);
    x_iter(:, iter+1) = x_iter(:, iter) + alpha_iter * p_iter(:, iter);
    r_iter(:, iter+1) = r_iter(:, iter) - alpha_iter * M * p_iter(:, iter);
    
    % Construct the next base vector p(:, iter+1)
    beta = norm(r_iter(:, iter+1))^2 / norm(r_iter(:, iter))^2;
    p_iter(:, iter+1) = r_iter(:, iter+1) + beta * p_iter(:, iter);
end

x_res = x_iter(:, iter);
% this is error 2
error = norm(r_iter(:, iter));

intermediates.iter_num = iter;
intermediates.x_iter = x_iter(:, 1:iter);
intermediates.r_iter = rn_iter(:, 1:iter)';
intermediates.p_iter = p_iter(:, 1:iter);
intermediates.err2_iter = err2_iter(:, 1:iter)';
intermediates.errA_iter = errA_iter(:, 1:iter)';

end

