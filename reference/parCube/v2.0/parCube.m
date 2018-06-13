%Vagelis Papalexakis, 2012
%School of Computer Science, Carnegie Mellon University
%Implementation of ParCube Non-negative PARAFAC decomposition for
%memory-resident tensors
function [A B C lambda] = parCube(X,F,sample_factor,times,nonneg)

if nargin == 4
    nonneg = 0;
end
mypath = pwd;


p = 0.55;
s = size(X); I = s(1); J = s(2); K = s(3);
A = sparse(I,F); B = sparse(J,F);C = sparse(K,F);
for i = 1:times
    if i == 1
        [Xs idx_i idx_j idx_k] = parCube_core(X,sample_factor);
        fixed_i = idx_i(1:ceil(p*length(idx_i)));
        fixed_j = idx_j(1:ceil(p*length(idx_j)));
        fixed_k = idx_k(1:ceil(p*length(idx_k)));
        fixed{1} = fixed_i; fixed{2} = fixed_j; fixed{3}=fixed_k;
    else 
        [Xs idx_i idx_j idx_k] = parCube_core(X,sample_factor,fixed);
    end
    
    if nnz(Xs)>0
        if nonneg
            factors = cp_nmu(Xs,F);
        else
            factors = cp_als(Xs,F);
        end
        As=sparse(factors.U{1});Bs=sparse(factors.U{2});Cs=sparse(factors.U{3});%lambda = factors.lambda;As = As*diag(lambda);
        lambda(:,i) = factors.lambda;
        As = As*diag(lambda(:,i).^(1/3)); Bs = Bs*diag(lambda(:,i).^(1/3)); Cs = Cs*diag(lambda(:,i).^(1/3));  
    else
       As = 0; Bs = 0; Cs = 0; 
    end
    
    
    Atemp = sparse(I,F); Btemp = sparse(J,F);Ctemp = sparse(K,F);
    Atemp(idx_i,:) = As; Btemp(idx_j,:) = Bs; Ctemp(idx_k,:) = Cs;
    
     %now, normalize the common part of every column
    for f = 1:F
       norm_a = norm(Atemp(fixed_i,f),2);Atemp(:,f) = Atemp(:,f)/norm_a; lambda(f,i) = norm_a;
       norm_b = norm(Btemp(fixed_j,f),2);Btemp(:,f) = Btemp(:,f)/norm_b; lambda(f,i) = lambda(f,i)*norm_b;
       norm_c = norm(Ctemp(fixed_k,f),2);Ctemp(:,f) = Ctemp(:,f)/norm_c; lambda(f,i) = lambda(f,i)*norm_c;
    end
    
    %Do the merge
    if i ==1
        A = Atemp; B = Btemp; C = Ctemp;
    else
        valA=zeros(F,1);valB=zeros(F,1);valC=zeros(F,1);
        for f1 = 1:F
            for f2 = 1:F
                valA(f2) = Atemp(fixed_i,f1)'* A(fixed_i,f2);
                valB(f2) = Btemp(fixed_j,f1)'* B(fixed_j,f2);
                valC(f2) = Ctemp(fixed_k,f1)'* C(fixed_k,f2);
            end
            [junk idx] = max(valA);
            mask2 = find(A(:,idx)==0);
            A(mask2,idx) = A(mask2,idx) + Atemp(mask2,f1);%Update ONLY the zero values
            
            [junk idx] = max(valB);
            mask2 = find(B(:,idx)==0);
            B(mask2,idx) = B(mask2,idx) + Btemp(mask2,f1);%Update ONLY the zero values
            
            [junk idx] = max(valC);
            mask2 = find(C(:,idx)==0);
            C(mask2,idx) = C(mask2,idx) + Ctemp(mask2,f1);%Update ONLY the zero values
        end
    
    end

end

lambda = mean(lambda,2);
A = sparse(A*diag(lambda));
A(A<10^-8) = 0;
B(B<10^-8) = 0;
C(C<10^-8) = 0;

