function [y,B] = mscaling(dmatrix,r)
% y=mscaling(dmatrix,r); r:dimensionality of the new space 
% by anderson  

 
[m,n] = size(dmatrix);
 
a=-(1/2)*dmatrix;

ac=mean(a); al=mean(a'); aa=mean(mean(a));

    
	for i=1:n;
              for j=1:n;
       b(i,j)=a(i,j)-ac(i)-al(j)+aa;
              end
        end


   [v,d]=eig(b); % eigenvector is a column in v

  evalues=diag(d); 
  [h,hh]=sort(evalues);	
  
   for i=1:r; 
       c(:,i)= v(:,hh(n+1-i))* ((evalues(hh(n+1-i)))^(1/2)); 
       end

   y=c;     

 B=b;




