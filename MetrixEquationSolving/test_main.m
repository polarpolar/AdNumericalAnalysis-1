clear all; close all; clc

%% Generate Random Metrix or Load Metrix
%genEqu();
load equation.mat

err_thr = 1.0 * 10^-15;
max_iternum = 100;

errB_dgra_spd_mean_cg = zeros(1, 5);
err_dgra_spd_mean_cg = zeros(1, 5);
errA_dgra_spd_mean_cg = zeros(1, 5);
errB_dgra_spd_mean_lanc = zeros(1, 5);
err_dgra_spd_mean_lanc = zeros(1, 5);
errA_dgra_spd_mean_lanc = zeros(1, 5);
errB_dgra_spd_mean_minr = zeros(1, 5);
err_dgra_spd_mean_minr = zeros(1, 5);
errA_dgra_spd_mean_minr = zeros(1, 5);

errB_dgra_cg = zeros(max_iternum, 5);
err_dgra_cg = zeros(max_iternum, 5);
errA_dgra_cg = zeros(max_iternum, 5);
errB_dgra_lanc = zeros(max_iternum, 5);
err_dgra_lanc = zeros(max_iternum, 5);
errA_dgra_lanc = zeros(max_iternum, 5);
errB_dgra_minr = zeros(max_iternum, 5);
err_dgra_minr = zeros(max_iternum, 5);
errA_dgra_minr = zeros(max_iternum, 5);

errB_dgra_spd_cg = zeros(max_iternum-1, 5);
err_dgra_spd_cg = zeros(max_iternum-1, 5);
errA_dgra_spd_cg = zeros(max_iternum-1, 5);
errB_dgra_spd_lanc = zeros(max_iternum-1, 5);
err_dgra_spd_lanc = zeros(max_iternum-1, 5);
errA_dgra_spd_lanc = zeros(max_iternum-1, 5);
errB_dgra_spd_minr = zeros(max_iternum-1, 5);
err_dgra_spd_minr = zeros(max_iternum-1, 5);
errA_dgra_spd_minr = zeros(max_iternum-1, 5);


for i = 1:5
    M = A(:, :, i, 1);
    %[ M, b ] = JacobiPrepro( M, b );
    
%% C-G iteration to solve A_ij * x = b
[ iter_num_cg, x_cg, err_cg, interm_cg ] = cgIter( M, b, err_thr, max_iternum );
% Backward error (2-norm(r)) degration
errB_dgra_spd_mean_cg(:, i) = mean( interm_cg.r_iter(2:iter_num_cg) ./...
    interm_cg.r_iter(1:iter_num_cg-1) );
errB_dgra_cg(1:iter_num_cg, i) = interm_cg.r_iter(1:iter_num_cg);
errB_dgra_spd_cg(1:iter_num_cg-1, i) = interm_cg.r_iter(2:iter_num_cg) ./...
    interm_cg.r_iter(1:iter_num_cg-1);
% Forward error (2-norm) degration
err_dgra_spd_mean_cg(:, i) = mean( interm_cg.err2_iter(2:iter_num_cg) ./...
    interm_cg.err2_iter(1:iter_num_cg-1) );
err_dgra_cg(1:iter_num_cg, i) = interm_cg.err2_iter(1:iter_num_cg);
err_dgra_spd_cg(1:iter_num_cg-1, i) = interm_cg.err2_iter(2:iter_num_cg) ./...
    interm_cg.r_iter(1:iter_num_cg-1);
% Forward error (A-norm)  degration
errA_dgra_spd_mean_cg(:, i) = mean( interm_cg.errA_iter(2:iter_num_cg) ./...
    interm_cg.errA_iter(1:iter_num_cg-1) );
errA_dgra_cg(1:iter_num_cg, i) = interm_cg.errA_iter(1:iter_num_cg);
errA_dgra_spd_cg(1:iter_num_cg-1, i) = interm_cg.errA_iter(2:iter_num_cg) ./...
    interm_cg.errA_iter(1:iter_num_cg-1);

    figure(1)
    subplot(2, 2, 1); semilogy(err_dgra_cg(:, i)); hold on
    if i == 5
        title('Forward error (deltaX) in C-G');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end
    subplot(2, 2, 2); plot(err_dgra_spd_cg(:, i)); hold on
    if i == 5
        title('Forward converge speed in C-G');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end
    subplot(2, 2, 3); semilogy(errB_dgra_cg(:, i)); hold on
    if i == 5
        title('Backward error (r) in C-G');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end
    subplot(2, 2, 4); plot(errB_dgra_spd_cg(:, i)); hold on
    if i == 5
        title('Backward converge speed in C-G');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end


%% Lanczos iteration to solve A_ij * x = b
[ iter_num_lanc, x_lanc, err_lanc, interm_lanc ] = lanczosIter( M, b, err_thr, max_iternum );
% Backward error (2-norm(r)) degration
errB_dgra_spd_mean_lanc(:, i) = mean( interm_lanc.r_iter(2:iter_num_lanc) ./...
    interm_lanc.r_iter(1:iter_num_lanc-1) );
errB_dgra_lanc(1:iter_num_lanc, i) = interm_lanc.r_iter(1:iter_num_lanc);
errB_dgra_spd_lanc(1:iter_num_lanc-1, i) = interm_lanc.r_iter(2:iter_num_lanc) ./...
    interm_lanc.r_iter(1:iter_num_lanc-1);
% Forward error (2-norm) degration
err_dgra_spd_mean_lanc(:, i) = mean( interm_lanc.err2_iter(2:iter_num_lanc) ./...
    interm_lanc.err2_iter(1:iter_num_lanc-1) );
err_dgra_lanc(1:iter_num_lanc, i) = interm_lanc.err2_iter(1:iter_num_lanc);
err_dgra_spd_lanc(1:iter_num_lanc-1, i) = interm_lanc.err2_iter(2:iter_num_lanc) ./...
    interm_lanc.r_iter(1:iter_num_lanc-1);
% Forward error (A-norm)  degration
errA_dgra_spd_mean_lanc(:, i) = mean( interm_lanc.errA_iter(2:iter_num_lanc) ./...
    interm_lanc.errA_iter(1:iter_num_lanc-1) );
errA_dgra_lanc(1:iter_num_lanc, i) = interm_lanc.errA_iter(1:iter_num_lanc);
errA_dgra_spd_lanc(1:iter_num_lanc-1, i) = interm_lanc.errA_iter(2:iter_num_lanc) ./...
    interm_lanc.errA_iter(1:iter_num_lanc-1);

    figure(2)
    subplot(2, 2, 1); semilogy(err_dgra_lanc(:, i)); hold on
    if i == 5
        title('Forward error (deltaX) of Lanczos');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end
    subplot(2, 2, 2); plot(err_dgra_spd_lanc(:, i)); hold on
    if i == 5
        title('Forward converge speed in Lanczos');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end
    subplot(2, 2, 3); semilogy(errB_dgra_lanc(:, i)); hold on
    if i == 5
        title('Backward error (r) in Lanczos');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end
    subplot(2, 2, 4); plot(errB_dgra_spd_lanc(:, i)); hold on
    if i == 5
        title('Backward converge speed in Lanczos');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end

%% MINRES iteration to solve A-ij * x = b
[ iter_num_minr, x_minr, err_minr, interm_minr ] = minresIter( M, b, err_thr, max_iternum );
% Backward error (2-norm(r)) degration
errB_dgra_spd_mean_minr(:, i) = mean( interm_minr.r_iter(2:iter_num_minr) ./...
    interm_minr.r_iter(1:iter_num_minr-1) );
errB_dgra_minr(1:iter_num_minr, i) = interm_minr.r_iter(1:iter_num_minr);
errB_dgra_spd_minr(1:iter_num_minr-1, i) = interm_minr.r_iter(2:iter_num_minr) ./...
    interm_minr.r_iter(1:iter_num_minr-1);
% Forward error (2-norm) degration
err_dgra_spd_mean_minr(:, i) = mean( interm_minr.err2_iter(2:iter_num_minr) ./...
    interm_minr.err2_iter(1:iter_num_minr-1) );
err_dgra_minr(1:iter_num_minr, i) = interm_minr.err2_iter(1:iter_num_minr);
err_dgra_spd_minr(1:iter_num_minr-1, i) = interm_minr.err2_iter(2:iter_num_minr) ./...
    interm_minr.r_iter(1:iter_num_minr-1);
% Forward error (A-norm) degration
errA_dgra_spd_mean_minr(:, i) = mean( interm_minr.errA_iter(2:iter_num_minr) ./...
    interm_lanc.errA_iter(1:iter_num_minr-1) );
errA_dgra_minr(1:iter_num_minr, i) = interm_minr.errA_iter(1:iter_num_minr);
errA_dgra_spd_minr(1:iter_num_minr-1, i) = interm_minr.errA_iter(2:iter_num_minr) ./...
    interm_minr.errA_iter(1:iter_num_minr-1);

    figure(3)
    subplot(2, 2, 1); semilogy(errA_dgra_minr(:, i)); hold on
    if i == 5
        title('Forward error (deltaX) in MINRES');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end
    subplot(2, 2, 2); plot(errA_dgra_spd_minr(:, i)); hold on
    if i == 5
        title('Forward converge speed in MINRES');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end
    subplot(2, 2, 3); semilogy(errB_dgra_minr(:, i)); hold on
    if i == 5
        title('Backward error (r) in MINRES');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end
    subplot(2, 2, 4); plot(errB_dgra_spd_minr(:, i)); hold on
    if i == 5
        title('Backward converge speed in MINRES');
        legend('D1', 'D2', 'D3', 'D4', 'D5');
    end

end


%iter_size = min([iter_num_cg, iter_num_lanc, iter_num_minr]);

% figure 
% plot(interm_cg.err2_iter(1:iter_size)); hold on
% plot(interm_lanc.err2_iter(1:iter_size)); hold on
% plot(interm_minr.err2_iter(1:iter_size)); hold on
% title('err2 of each iteration');
% legend('c-g', 'lanczos', 'minres');
% 
% figure 
% plot(interm_cg.errA_iter(1:iter_size)); hold on
% plot(interm_lanc.errA_iter(1:iter_size)); hold on
% plot(interm_minr.errA_iter(1:iter_size)); hold on
% title('errA of each iteration');
% legend('c-g', 'lanczos', 'minres');


