clear all; close all; clc;
figure(1)
grid on;
n=1000;
%% Def
lamda=sort(rand(n,1),'descend');
D=diag(lamda);Q=orth(rand(n,n));b=rand(n,1);A=Q*D*Q';
P=zeros(n,n);R=zeros(n,n);X=zeros(n,n);alpha=zeros(n,1);beta=alpha;
%% CG
r_cg=CG(n,A,b);
semilogy(r_cg);
grid on
hold on
%% Lanczos
loss_lans=Lanczos(n,A,b);
semilogy(loss_lans);
xlabel('µü´ú´ÎÊý k');ylabel('||r_k||');
legend('CG','Lanczos');

