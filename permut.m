function   y=permut(x)
% function   y=permut(x)
%  random permutation of the elements of x vector


r=randn(1,x);
[sr,isr]=sort(r);

y=isr;

