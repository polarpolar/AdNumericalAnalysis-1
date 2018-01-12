clear all; close all; clc;
n=1000;
for m=[0,10,100,500]
    figure()
    grid on;
    %% Def
    lamda(1:m)=-rand(m,1);lamda(m+1:n)=rand(n-m,1);
    lamda=sort(lamda,'descend');
    D=diag(lamda);Q=orth(rand(n,n));A=Q*D*Q';b=rand(n,1);
    P=zeros(n,n);R=zeros(n,n);X=zeros(n,n);alpha=zeros(n,1);beta=alpha;
    %% Lanczos
    loss_lans=Lanczos_n(n,A,b);
    semilogy(loss_lans);
    hold on;
    %% MINRES
    loss_lans=MINRES(n,A,b);
    semilogy(loss_lans);
    hold on;
    title(['负特征值数目m=',num2str(m)]);xlabel('迭代次数 k');ylabel('||r_k||');
    legend('Lanczos','MINRES');
end
