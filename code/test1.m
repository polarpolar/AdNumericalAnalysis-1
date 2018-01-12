clear all; close all; clc;
figure(1)
grid on;
n=1000;
for k=[10,100,1000];
    %% Def
    lamda=linspace(1,k,n);
    D=diag(lamda);Q=orth(rand(n,n));b=rand(n,1);A=Q*D*Q';
    P=zeros(n,n);R=zeros(n,n);X=zeros(n,n);alpha=zeros(n,1);beta=alpha;
    %% CG
    r_cg=CG(n,A,b);
    semilogy(r_cg);
    hold on
end
xlabel('µü´ú´ÎÊý k');ylabel('||r_k||');
legend('k=10','k=100','n=1000');
grid on;
