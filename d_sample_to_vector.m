function  y=d_sample_to_vector(X,Y)

% distance=d_sample_to_vector(sample,vector)
% distances (squarred euclidean to a reference-vector)

[m,n]=size(X);


y= diag([X-ones(m,1)*Y]*[X-ones(m,1)*Y]');