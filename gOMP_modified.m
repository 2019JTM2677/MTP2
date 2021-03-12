function [x]=gOMP_modified(K,y,A,N)
%A: Sensing matrix
%y: measurement vector
%K: Sparsity
%N: No. of indices chosen in each iteration
[m,n] = size (A) ;
%x = zeros (n,1) ;
if (N > K) || (N > m/K)
 fprintf("Give correct value for N");
 return 
end

k=0;
Residue =[];
B=[];      
lambda={};% support set
e = 10^-5;
N_d=N+1;
L=nchoosek(N_d,N);

while k< min(K,m/N)
    k=k+1;
    index=[];d=[]; indices={};
    if k==1
        [~ ,index] = maxk(abs(A' * y),N_d);
        d = nchoosek(index,N);
        for j=1:L
            indices{end+1}=d(j,:);
            %display(indices{end})
        end
    else
       
        for i=1:size(lambda,2)
            if norm(Residue(:,i)) > e
                [~ ,index(:,i)] = maxk(abs(A' * Residue(:,i)),N_d);
                %display(abs(A' * Residue(:,i)));
                %display(val)
                %display(index(:,i))
                d=[repmat(lambda{i},L,1) nchoosek(index(:,i),N)];
                for j=1:L
                    indices{end+1}=d(j,:);
                    %display(indices{end})
                end
            else 
                d = lambda{i};
                indices{end+1}=d;
                %display(indices{end})
            end
            %display(d); display(size(indices));
            
        end
    end
    
    lambda = indices;
    %lambda = [lambda index];
    %display(lambda);
    B=[];
    x_cap=[];x_cap_save={};
    for i=1:size(lambda,2)
        
        B= A(:,lambda{i});
        x_cap = pinv(B)*y;
        %x_cap = B_mod\y_mod;
        [warnMsg, warnID] = lastwarn();
        if(isempty(warnID))
        else
            display(lambda{i})
            display(y)
            display(A)
        end
        x_cap_save{i}=x_cap;
        
        
        Residue(:,i) = y - B*x_cap;
        %display(Residue)
    end
end
x = zeros (n,size(lambda,2)) ;

for i=1:size(lambda,2)
    
    x(lambda{i},i)= x_cap_save{i};
end
end