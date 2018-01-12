clear all; close all; clc

% load data
load result.mat;

% select five metrics with each method to draw the figures.
selection = [...
    1, 1;
    1, 2;
    1, 3;
    1, 4;
    1, 5;
    ];

for alg = 1 : 3
    switch alg
        case 1
            name = 'C-G';
        case 2
            name = 'Lanczos';
        case 3
            name = 'MINRES';
    end
    figure;
    for sel = 1 : 5
        errA_dgra_tmp(1, :) = errA_iter(selection(sel, 1), selection(sel, 2),...
            alg, 2:end);
        plot(errA_dgra_tmp); hold on;
    end
    title(['The convergency curve of ', name]);
    legend('M1', 'M2', 'M3', 'M4', 'M5');
end

clear errA_dgra_tmp;
clear name;
clear sel;
clear alg;
%close all;
