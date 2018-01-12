function [ iter, x_res, error, intermediates ] = minresIter( M, b, err_thr, max_iter )
%MINRESITER this is the MINRES solution of equation solving.
%   The MINRES solution of metrix equation solving. 

% check arguments
[Msize, ~] = size(M);
if nargin <= 3; max_iter = Msize + 1; end

%% MINRES Method to solve Equation
% accurate result
x_acc = M \ b;

% variables
x_iter = zeros(Msize, max_iter);
y_iter = zeros(Msize, max_iter);
r_iter = zeros(1, max_iter);
err2_iter = zeros(1, max_iter);
errA_iter = zeros(1, max_iter);

% initial value
x_init = zeros(Msize, 1);
beta_init = norm(b - M * x_init);

%   Give the complete lanczos process. All the intermediate metrics are
%   defined as the k-th sub-metrix of Q.
[ Q, T ] = lanczosProc( M, b, Msize );
for iter = 1:max_iter
    % Find the base metrix of the k-th sub-space. The base metrix is
    % constructed by the k-th sub-metrix of lanczos metrix Q, and T. 
    % get Q^(k+1) & T^(k) with Lanczos Process
	Q_iter = Q(:, 1:iter);
    if iter == Msize
        T_iter = [T(1:iter, 1:iter); zeros(1, iter)];
    else
        T_iter = T(1:iter+1, 1:iter);
    end
    
    % QR decompose T with GIVENs.
    %       P is (k+1) * (k+1), and R is (k+1) * k.
    [ P, R ] = qr(T_iter); %givensCvrt( T, k );
    % If R is singular, continue ... (Actually, it shouldn't happen.)
    if rank(R) < iter
        disp('MINRES: QR decomposion is singular.');
        continue;
    end
    
    %  solve simplified equation
    %       R^ * x = b;
    %       R^ = R(1:k, :)
    %       beta * P' * e1 = [b', b(k+1)]'
    %       b(k+1) is error-r_k
    e1 = zeros(iter+1, 1);     e1(1) = beta_init;
    b = P' * e1;    
    y_iter(1:iter, iter)  = R(1:iter, :) \ b(1:iter);
	x_iter(:, iter) = Q_iter(:, 1:iter) * y_iter(1:iter, iter)  + x_init;

    r_iter(iter) = abs(b(iter+1));
    
    % Forward Error
    err2_iter(iter) = norm( x_iter(:, iter) - x_acc );
    errA_iter(iter) = Anorm( x_iter(:, iter) - x_acc, M );
    
    % Termination Judgement
    if r_iter(iter) < err_thr
        break
    end  
end

x_res = x_iter(:, iter);
error = r_iter(iter);

intermediates.iter_num = iter;
intermediates.x_iter = x_iter(:, 1:iter);
intermediates.r_iter = r_iter(1:iter)';
intermediates.err2_iter = err2_iter(1:iter)';
intermediates.errA_iter = errA_iter(1:iter)';

end

