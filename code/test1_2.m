clear all; close all; clc;
figure(1)
grid on;
n=1000;lamda=linspace(1,n,n)/n;
for k=[1,100,10000];
    lamda_=lamda;
    %lamda_(1)=1/n/k;
    lamda_(n)=1*k;
    %% Def
    D=diag(lamda_);Q=orth(rand(n,n));b=rand(n,1);A=Q*D*Q';
    P=zeros(n,n);R=zeros(n,n);X=zeros(n,n);alpha=zeros(n,1);beta=alpha;
    %% CG
    r_cg=CG(n,A,b);
    semilogy(r_cg);
    hold on
end
xlabel('迭代次数 k');ylabel('||r_k||');
legend('特征值正态分布','最大、最小特征值扩大缩小10^2倍','最大、最小特征值扩大缩小10^4倍');
grid on





