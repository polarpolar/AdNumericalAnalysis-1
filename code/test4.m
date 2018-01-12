clear all; close all; clc;
figure(1)
grid on;
n=1000;
for m=[10,100,500,1000]
    %% Def
    lamda=sort(rand(n,1),'descend');
    D=diag(lamda);Q=orth(rand(n,n));A=Q*D*Q';
    seq=randperm(n);a=rand(m,1);b=Q(:,seq(1:m))*a/norm(a);
    P=zeros(n,n);R=zeros(n,n);X=zeros(n,n);alpha=zeros(n,1);beta=alpha;
    %% Lanczos
    loss_lans=Lanczos(n,A,b);
    semilogy(loss_lans);
    hold on;
end
grid on;
xlabel('µü´ú´ÎÊý k');ylabel('||r_k||');
legend('b=10','b=100','b=500','b=1000');