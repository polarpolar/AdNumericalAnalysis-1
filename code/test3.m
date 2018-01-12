clear all; close all; clc;
figure(1)
grid on;
n=1000;
for m=[10,100,500,990,1000]
    %% Def
    lamda(1:m)=rand(m,1);
    if(m<1000)
        lamda(m+1:n)=lamda(ceil(m*rand(1,n-m)));
    end
    lamda=sort(lamda,'descend');
    D=diag(lamda);Q=orth(rand(n,n));b=rand(n,1);A=Q*D*Q';
    P=zeros(n,n);R=zeros(n,n);X=zeros(n,n);alpha=zeros(n,1);beta=alpha;
    %% Lanczos
    loss_lans=Lanczos(n,A,b);
    semilogy(loss_lans);
    hold on;
end
grid on;
xlabel('µü´ú´ÎÊý k');ylabel('||r_k||');
legend('m=10','m=100','m=500','m=990','m=1000');