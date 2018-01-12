function [ Q, R ] = givensCvrt( T, k )
%givensCvrt Summary of this function goes here
%   Detailed explanation goes here

T_iter = zeros(k, k-1, k-1);
T_iter(:, :, 1) = T(:, 1:k-1);    % T: k * k-1

G = zeros(k, k, k);
Q = diag(ones(1,k));
for i = 1:k-1
    G(:, :, i) = diag(ones(1, k));
    tmp = sqrt( T_iter(i, i, i).^2 + T_iter(i+1, i, i).^2 );
    G(i, i, i) = T_iter(i, i, i) / tmp;
    G(i+1, i+1, i) = T_iter(i, i, i) / tmp;
    G(i, i+1, i) = T_iter(i+1, i, i) / tmp;
    G(i+1, i, i) = -T_iter(i+1, i, i) / tmp;
    T_iter(:, :, i+1) = G(:, :, i) * T_iter(:, :, i);
    Q = G(:, :, i) * Q;
end

Q = Q';
R = T_iter(:, :, i+1);
%R = R .* (R > 0.001);

end

