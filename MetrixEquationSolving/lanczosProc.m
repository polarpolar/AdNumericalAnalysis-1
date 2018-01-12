function [ Q, T ] = lanczosProc( M, b, Msize )
%LANCZOSPROC is the lanczos process to generate the target metrix T and the
%   transformer base metrix Q based on M. 
%   In lanczos process, Q is sole only if the first vecter is the same.
%   Target on solving equations, we set the first vector (base) base on the
%   equation of Mx = b.

%% Lanczos Process ( Q'AQ = T, find T)
x_lanc = zeros(Msize, Msize+1);
r_lanc = zeros(Msize, Msize+1);
v_lanc = zeros(Msize, Msize+1);
alpha_lanc = zeros(1, Msize+1);
beta_lanc = zeros(1, Msize+1);

% initial value
x_lanc(:, 1) = zeros(Msize, 1);
r_lanc(:, 1) = b - M * x_lanc(:, 1);
v_lanc(:, 2) = r_lanc(:, 1) / sqrt(sum( r_lanc(:,1).^2 ));
beta_lanc(1) = 0;

% iteration
for k = 2 : Msize + 1
    % find the length of the k-th base vector in M * v_lanc(:, k)
    alpha_lanc(:, k) = v_lanc(:, k)' * M * v_lanc(:, k);
    
    % remove the weights of <=k-th base vectors, and the remained is
    % the (k+1)-th base vector.
    f = M * v_lanc(:, k) - beta_lanc(:, k-1) * v_lanc(:, k-1)...
        - alpha_lanc(:, k) * v_lanc(:, k);
    
    % Even if the space is not complete, we do not need more dimensions.
    if k == Msize + 1; break; end
    
    % check if the space is complete
    if norm( f )
        % if the remained is not zero, the (k+1)-th base vector exists.
        % Get the weight
        beta_lanc(:, k) = norm( f );
        % Get the direction (base vector)
        v_lanc(:, k+1) = f/beta_lanc(:, k);
    else
        % if the remained is zero, then we find the complete space bases.
        %fprintf('Benign quit @ %d!\n', k);
        break
    end
    
end

% alpha: the length is Msize
alpha_lanc = alpha_lanc(2:Msize+1);
% beta: the length is Msize - 1
beta_lanc = beta_lanc(2:Msize);

Q = v_lanc(:, 2:Msize+1);
% T = (Q' * M * Q);% .* (Q' * M * Q > 0.0001);
T = diag(alpha_lanc) + diag(beta_lanc, 1) + diag(beta_lanc, -1);

end

