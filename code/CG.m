function [ r ] = CG( n,A,b )
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
P=zeros(n,n);R=zeros(n,n);X=zeros(n,n);alpha=zeros(n,1);beta=alpha;
%% CG
% initial
x0=zeros(n,1);
r0=b-A*x0;
P(:,1)=r0;
alpha(1)=(r0'*r0)/(P(:,1)'*A*P(:,1));
X(:,1)=x0+alpha(1)*P(:,1);
R(:,1)=b-A*X(:,1);
beta(2)=(R(:,1)'*R(:,1))/(r0'*r0);
P(:,2)=R(:,1)+beta(2)*P(:,1);
%iter
R_1=(R(:,1)'*R(:,1));AP=A*P(:,2);R_2=R_1;
for i=2:n
    alpha(i)=R_2/(P(:,i)'*AP);
    X(:,i)=X(:,i-1)+alpha(i)*P(:,i);
    R(:,i)=R(:,i-1)-alpha(i)*AP;R_1=(R(:,i)'*R(:,i));
    beta(i+1)=R_1/R_2;
    P(:,i+1)=R(:,i)+beta(i+1)*P(:,i);
    AP=A*P(:,i+1);R_2=R_1;
end
r=(sqrt(sum(R.*R)));
end

