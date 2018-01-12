function [ loss ] = Lanczos_n( n,A,b)
alpha=zeros(n,1);beta=zeros(n,1);Q=zeros(n,n);R=zeros(n,n);
T=zeros(n,n);loss=zeros(n,1);
% initial
x0=zeros(n,1);
r0=b-A*x0;
nr=norm(r0);
Q(:,1)=r0/nr;

r=A*Q(:,1);
alpha(1)=Q(:,1)'*r;
r=r-alpha(1)*Q(:,1);
beta(1)=norm(r);
Q(:,2)=r/beta(1); 
T(1,1)=alpha(1);T(1,2)=beta(1);T(2,1)=beta(1);
y=nr/T(1,1);
loss(1)=beta(1)*y(1);
% iter
for j=2:n
    r=A*Q(:,j)-beta(j-1)*Q(:,j-1);
    alpha(j)=Q(:,j)'*r;
    r=r-alpha(j)*Q(:,j);
    beta(j)=norm(r);
    if (beta(j)==0)
    begin
        break;
    end
    Q(:,j+1)=r/beta(j); 
    T(j,j)=alpha(j);T(j,j+1)=beta(j);T(j+1,j)=beta(j);
    y=(T(1:j,1:j))\[1;zeros(j-1,1)];
    loss(j)=beta(j)*y(j);
end
loss=abs(loss);
end