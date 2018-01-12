function [ y ] = S2( L,b)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
n=size(L,1);
y=zeros(n,1);
y(n)=b(n)/L(n,n);
for i=n-1:1
   y(i)=(b(i)-L(i,i+1:n)*y(i+1:n))/L(i,i); 
end
end

