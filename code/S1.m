function [ y ] = S1( L,nr )
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
n=size(L,1);
y=zeros(n,1);
y(1)=nr/L(1,1);
for i=2:n
   y(i)=(-L(i,1:i-1)*y(1:i-1))/L(i,i); 
end
end

