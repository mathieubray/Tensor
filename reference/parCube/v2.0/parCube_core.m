%Vagelis Papalexakis, 2012
%School of Computer Science, Carnegie Mellon University
%Core sampling function for the ParCube algorithm, for memory resident
%tensors.
function [Xs idx_i idx_j idx_k Xma Xmb Xmc] = parCube_core(X,sample_factor,fixed_set)

if numel(sample_factor)>1
    s1 = sample_factor(1); 
    s2 = sample_factor(2);
    s3 = sample_factor(3);
else
    s1 = sample_factor; 
    s2 = sample_factor;
    s3 = sample_factor;
end

if nargin == 3
   fixed_i = fixed_set{1};
   fixed_j = fixed_set{2};
   fixed_k = fixed_set{3};
else
    fixed_i = [];
    fixed_j = [];
    fixed_k = [];
end
    
negative_values = find(X<0);
if ~isempty(negative_values)
    X_abs = double(X);
    X_abs = abs(X_abs);
    X_abs = sptensor(X_abs);
else
   X_abs = X; 
end
s = size(X); I = s(1); J = s(2); K = s(3);


Xma = collapse(X_abs,[2 3]);
Xmb = collapse(X_abs,[1 3]);
Xmc = collapse(X_abs,[1 2]);


idx_i = find(Xma~=0);
idx_i = randsample_weighted(idx_i,I/s1,full(Xma(idx_i)));
if ~isempty(fixed_i)
    idx_i = union(idx_i,fixed_i);
end

idx_j = find(Xmb~=0);
idx_j = randsample_weighted(idx_j,J/s2,full(Xmb(idx_j)));
if ~isempty(fixed_j)
    idx_j = union(idx_j,fixed_j);
end

idx_k = find(Xmc~=0);
idx_k = randsample_weighted(idx_k,K/s3,full(Xmc(idx_k)));
if ~isempty(fixed_k)
    idx_k = union(idx_k,fixed_k);
end

Xs = X(idx_i,idx_j,idx_k);