% Generalized OMP
% https://arxiv.org/pdf/1111.6664.pdf

function [x]=gOMP(K,y,A,N)
%A: Sensing matrix
%y: measurement vector
%K: Sparsity
%N: No. of indices chosen in each iteration
[m,n] = size (A) ;
x = zeros (n,1) ;
if (N > K) || (N > m/K)
 fprintf("Give correct value for N");
 return 
end

k=0;
Residue =y;
B=[];      
lambda=[];% support set
e = 10^-5;

while k< min(K,m/N) && norm(Residue) > e
    k=k+1;
    [~ ,index] = maxk(abs(A' * Residue),N);
    %display(abs(A' * Residue))
    %display(val)
    %display(index)
    lambda = union(lambda,index);
    %display(lambda);
    B= A(:,lambda);
    %display(B);
    B_mod = B'*B;
    y_mod = B'*y;
   % x_cap = B_mod\y_mod;
    x_cap = pinv(B)*y;
    %display(x_cap)
    Residue = y - B*x_cap;
    %display(Residue)
end
x(lambda)=x_cap;
end