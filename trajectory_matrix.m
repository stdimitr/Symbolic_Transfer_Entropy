function    U = trajectory_matrix(x,l)

% U = trajectory_matrix(x,l)
%        
%    a signal-vector [x(1), x(2),.....x(T)]
%                                            
%    is converted to an array: U=[x(1),x(2),...x(l); 
%                               x(2),x(3),...x(l+1);                    
%                               ..............,,;                                                                  
%                               .............x(T)]                
%
%   l is the length of the segments          
%



T=length(x);
       
U=[];
    
    for i=l:T ;
      U=[U;  x(i-l+1 :i )];
    end  
   