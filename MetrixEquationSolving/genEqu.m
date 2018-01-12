function genEqu(  )
%GENEQU Summary of this function goes here
%   Detailed explanation goes here

%% Generate Source Matrix
% provide eigenvalue (100*100 diagonal matrix)
D = zeros(100, 100, 5);
% for i = 1:5
%     D(:,:,i) = 100*diag(rand(1,100), 0);
% end
D(:,:,1) = 10 * diag( sort([10, rand(1,98)*0.9+0.1, 0.1], 'descend') );
D(:,:,2) = 10 * diag( sort([10, rand(1,98)/2+9.5, 0.1], 'descend') );
D(:,:,3) = 10 * diag( sort([10, rand(1,98)*9+1, 0.1], 'descend') );
D(:,:,4) = 10 * diag( sort([10, (rand(1,98)*0.9+0.1)*10, 0.1], 'descend') );
D(:,:,5) = 10 * diag( sort([10, rand(1,98)+9, 9], 'descend') );

% generate 10 random 100*100 matrices
M = 100*rand(100, 100, 10);
Q = zeros(100, 100, 10);
for j = 1:10
    [Q(:, :, j), ~] = qr(M(:,:, j));
end
% generate 50 origin matrices
A = zeros(100, 100, 5, 10);
condNum = zeros(5, 10);    % condition number of A
condRatio = zeros(5, 10);
for i = 1:5
    for j = 1:10
        A(:, :, i, j) = Q(:, :, j)' * D(:, :, i) * Q(:, :, j);
        condNum(i, j) = cond(A(:, :, i, j), 2);
        condRatio(i, j) = (condNum(i, j)-1)/(condNum(i, j)+1);
    end
end

%% Solve Equation Ax = b
% provide a 'b'
b = rand(100, 1);

M = A(:, :, 1, 1);

save('equation.mat', 'M', 'b', 'A', 'D', 'Q');

end

