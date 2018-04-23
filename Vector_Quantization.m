function [prototypes,class_indicators,Average_Error,Convergence_Index]=Vector_Quantization(X,k,iterations)

% [prototypes,class_indicators,Average_Error,Convergence_Index]=Vector_Quantization(X,k,iteration_factor)
% 
%  Neural-Gas Vector Quantization Algorithm 
%  X is the input set of patterns, k the size of the code-book (i.e. number of prototypes/centroids)
%  the iteration_factor controls the number of iterations : == (iteration_factor) x (size of the input sample)
%
%  prototypes: tabulates the k code-vectors
%  class_indicators: tabulates the labels that assign each vector Xi to the nearest prototype
%  Average_Error is an index of performance: it is the average Distortion induced by the adopted coding scheme
%  Convergence_Index indicates the improvement,with respect to the initial/random selection of prototypes ,
%  achieved with the iterative-execution of the basic adaptation-step 

[N,p]=size(X);


%%%%% initialization %%%%%%%%%%%

rindex=permut(N);

fl=floor(N/k); 

      if fl>=2
         for i=1:k
         rr=rindex((i-1)*fl+1:(i)*fl);
         prot(i,:)=mean(X(rr,:));
         end

      else
         prot=X(rindex(1:k),:);
      end

%%%%%%%% initial coding error %%%%%%%

for i=1:N
d=d_sample_to_vector(prot,X(i,:));
[error(i)]=min(d);
end
initial_Average_Error=mean(error);
%______________________________________





tmax=iterations*N; 
rr=[]; 
for i=1:iterations, 
    rr=[rr;permut(N)']; 
end

li=0.3*k; lf=0.01; ei=0.5; ef=0.005; 
lt=li; et=ei;



for i=1:tmax;
    u=X(rr(i),:);
    du=d_sample_to_vector( prot,u);
   [sdu,ordering_list]=sort(du);
   [ignore,order]=sort(ordering_list); order=order-1;
    hl=exp(-order/lt);


   % for i=1:k; dprot(i,:)=et * hl(i) * (u-prot(i,:)); end
   
   dprot= et*repmat(hl,1,p).* (repmat(u,k,1)-prot);
   prot=prot+dprot;

   lt=li*(lf/li)^(i/tmax); et=ei*(ef/ei)^(i/tmax);
end

prototypes=prot;

for i=1:N
    d=d_sample_to_vector(prototypes,X(i,:));
    [error(i),indicator(i)]=min(d);
end

class_indicators=indicator;
Average_Error=mean(error);
Convergence_Index = abs(Average_Error-initial_Average_Error)/initial_Average_Error;




