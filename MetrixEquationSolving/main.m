clear all; close all; clc

%% Generate Random Metrix or Load Metrix
%genEqu();
load equation.mat

err_thr = 1.0 * 10^-20;
max_iternum = 100;
err2_final = zeros(5, 10, 3);
errA_dgra = zeros(5, 10, 3);
errA_iter = zeros(5, 10, 3, max_iternum);
condD = zeros(5, 1);

% select Aij, i = 1, ..., 5, j = 1, ..., 10. 'i' is the index of D, and j
% is the index of Q.
for i = 1:5
    % condition number
    condD(i) = cond(D(:,:,i));
    % judge if Metrixs are SPD
    isPD = sum ( diag( D(:,:,i) ) < zeros(100, 1) ) == 0;
    fprintf('-----------------------------------------------\n');
    if isPD
        fprintf('The %d-th diagnose: PD\n', i);
    else
        fprintf('The %d-th diagnose: un-PD\n', i);
    end
    for j = 1:10
        % the metrix is constructed by M = Qj' * Di * Qj
        fprintf('*******The %d-th Metrix*******\n', j);
        M = A(:, :, i, j);
        if isPD
            fprintf('---- C-G Iterating ----\n');
            % C-G iteration to solve A_ij * x = b
            [ iter_num_cg, x_cg, err_cg, interm_cg ] = cgIter( M, b, err_thr, max_iternum );
            % error degration
            errA_dgra_cg = mean( interm_cg.errA_iter(2:iter_num_cg) ./...
                interm_cg.errA_iter(1:iter_num_cg-1) );
            err2_final(i, j, 1)  = err_cg;
            errA_dgra(i, j, 1)  = errA_dgra_cg;
            errA_iter(i, j, 1, 2:iter_num_cg)  = interm_cg.errA_iter(2:iter_num_cg) ./...
                interm_cg.errA_iter(1:iter_num_cg-1);
            fprintf('---- C-G Completed! ----\n');
        end
        
        % Lanczos iteration to solve A_ij * x = b
        fprintf('---- Lanczos Iterating ----\n');
        [ iter_num_lanc, x_lanc, err_lanc, interm_lanc ] = lanczosIter( M, b, err_thr, max_iternum );
        % error degration
        errA_dgra_lanc = mean( interm_lanc.errA_iter(2:iter_num_lanc) ./...
            interm_lanc.errA_iter(1:iter_num_lanc-1) );
        err2_final(i, j, 2)  = err_lanc;
        errA_dgra(i, j, 2)  = errA_dgra_lanc;
        errA_iter(i, j, 2, 2:iter_num_lanc)  = interm_lanc.errA_iter(2:iter_num_lanc) ./...
            interm_lanc.errA_iter(1:iter_num_lanc-1);
        fprintf('---- Lanczos Completed ----\n');

        % MINRES iteration to solve A-ij * x = b
        fprintf('---- MINRES Iterating ----\n');
        [ iter_num_minr, x_minr, err_minr, interm_minr ] = minresIter( M, b, err_thr, max_iternum );
        % error degration
        errA_dgra_minr = mean( interm_minr.errA_iter(2:iter_num_minr) ./...
            interm_minr.errA_iter(1:iter_num_minr-1) );
        err2_final(i, j, 3)  = err_minr;
        errA_dgra(i, j, 3)  = errA_dgra_minr;
        errA_iter(i, j, 3, 2:iter_num_minr)  = interm_minr.errA_iter(2:iter_num_minr) ./...
            interm_minr.errA_iter(1:iter_num_minr-1);
        fprintf('---- MINRES Completed ----\n');

    end
end

% calculate the condition number of Ds


save('result.mat', 'A', 'b', 'D', 'Q', 'err2_final', 'errA_dgra', 'errA_iter', 'condD');


