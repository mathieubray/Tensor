function s = randsample_weighted(p,n,w)
%Vagelis Papalexakis, 2012
%School of Computer Science, Carnegie Mellon University
n = ceil(n);
nnz_w = nnz(w);
n = min([n nnz_w]);%this guard makes sure that we always have data to sample from
s = zeros(n,1);

for i = 1:n
    
if(strfind(version('-release'),'2009'))
    %http://www.mathworks.com/matlabcentral/newsreader/view_thread/255206
    if size(w,1) == 1
        w = w';
    end
    pp = w/sum(w);
    temp = [0;cumsum(pp)];
    e = min(temp,1);
    e(end) = 1;
    pp = diff(e);
    w = pp;
end
    s(i) = randsample(p,1,true,w);
    w(p==s(i))=0;
end