function [A B C D lambda] = file_parCube_4mode(FILENAME,dims,F,sample_factor,times)
%Vagelis Papalexakis, 2012
%School of Computer Science, Carnegie Mellon University
try
    rmdir('temp','s')
catch
    disp('Temporary storage is already formatted');
end

%END OF ARGUMENT STUFF    
if numel(sample_factor)>1
    s1 = sample_factor(1); 
    s2 = sample_factor(2);
    s3 = sample_factor(3);
    s4 = sample_factor(4);
else
    s1 = sample_factor; 
    s2 = sample_factor;
    s3 = sample_factor;
    s4 = sample_factor;
end
I = dims(1); J = dims(2); K = dims(3); L = dims(4);


Xma = sparse(I,1); Xmb = sparse(J,1); Xmc = sparse(K,1);Xmd = sparse(L,1);
MAXMEM = 2048;
% MAXMEM = 8192;

% tic
disp('calling external jar file')
system(sprintf('java -Xmx%dm -jar tensor-sampling-4mode.jar %s %d %d %d %d',MAXMEM,FILENAME,I,J,K,L));

FILENAME_IDX1 = sprintf('%s1.txt',FILENAME);
FILENAME_IDX2 = sprintf('%s2.txt',FILENAME);
FILENAME_IDX3 = sprintf('%s3.txt',FILENAME);
FILENAME_IDX4 = sprintf('%s4.txt',FILENAME);


% system('mkdir temp')
mkdir('temp');

disp('Created the marginals!')



A = sparse(I,F); B = sparse(J,F);C = sparse(K,F);D = sparse(L,F);
lambda = zeros(F,times);

p = 0.35;
%initialize the fixed set of indices
system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d 0',MAXMEM,FILENAME_IDX1,ceil(I/s1)));
idx_i = load(sprintf('%s.sample-0',FILENAME_IDX1));


system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d 0 ',MAXMEM,FILENAME_IDX2,ceil(J/s2)));
idx_j = load(sprintf('%s.sample-0',FILENAME_IDX2));

system(sprintf('java  -Xmx2048m -jar randsample-weighted.jar %s %d 0',FILENAME_IDX3,ceil(K/s3)));
idx_k = load(sprintf('%s.sample-0',FILENAME_IDX3));

system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d 0',MAXMEM,FILENAME_IDX4,ceil(L/s4)));
idx_l = load(sprintf('%s.sample-0',FILENAME_IDX4));

fixed_i = idx_i(1:ceil(p*length(idx_i)));
fixed_j = idx_j(1:ceil(p*length(idx_j)));
fixed_k = idx_k(1:ceil(p*length(idx_k)));
fixed_l = idx_l(1:ceil(p*length(idx_l)));


parfor t = 1:times
% for t = 1:times

disp(sprintf('Pass #%d',t))
    %now make one pass to the file and select only specific indices
    
    system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d %d',MAXMEM,FILENAME_IDX1,ceil(I/s1),t));
    idx_i = load(sprintf('%s.sample-%d',FILENAME_IDX1,t));
    system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d %d',MAXMEM,FILENAME_IDX2,ceil(J/s2),t));
    idx_j = load(sprintf('%s.sample-%d',FILENAME_IDX2,t));
    system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d %d',MAXMEM,FILENAME_IDX3,ceil(K/s3),t));
    idx_k = load(sprintf('%s.sample-%d',FILENAME_IDX3,t));
    system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d %d',MAXMEM,FILENAME_IDX4,ceil(L/s4),t));
    idx_l = load(sprintf('%s.sample-%d',FILENAME_IDX4,t));
    
    idx_i = union(idx_i,fixed_i);dlmwrite(sprintf('%s.sample-%d',FILENAME_IDX1,t),idx_i,'delimiter','','precision','%.0f');
    idx_j = union(idx_j,fixed_j);dlmwrite(sprintf('%s.sample-%d',FILENAME_IDX2,t),idx_j,'delimiter','','precision','%.0f');
    idx_k = union(idx_k,fixed_k);dlmwrite(sprintf('%s.sample-%d',FILENAME_IDX3,t),idx_k,'delimiter','','precision','%.0f');
    idx_l = union(idx_l,fixed_l);dlmwrite(sprintf('%s.sample-%d',FILENAME_IDX4,t),idx_k,'delimiter','','precision','%.0f');
    
    system(sprintf('java -Xmx%dm -jar create-sampled-tensor-4mode.jar %s.sample-%d %s.sample-%d %s.sample-%d %s.sample-%d %s %d',MAXMEM,FILENAME_IDX1,t,FILENAME_IDX2,t,FILENAME_IDX3,t,FILENAME_IDX4,t,FILENAME,t));
    dat = dlmread(sprintf('%s-%d',FILENAME,t),'\t');
    
    [junk tmp_i] = ismember(dat(:,1),idx_i);
    [junk tmp_j] = ismember(dat(:,2),idx_j);
    [junk tmp_k] = ismember(dat(:,3),idx_k);
    [junk tmp_l] = ismember(dat(:,4),idx_l);
    
    Xs = sptensor([tmp_i tmp_j tmp_k tmp_l],dat(:,end),[length(idx_i) length(idx_j) length(idx_k) length(idx_l)]);   
    
       factors = cp_nmu(Xs,F);
    As=sparse(factors.U{1});Bs=sparse(factors.U{2});Cs=sparse(factors.U{3});Ds=sparse(factors.U{4});%lambda = factors.lambda;As = As*diag(lambda);
    lambda(:,t) = factors.lambda;
    As = As*diag(lambda(:,t).^(1/4)); Bs = Bs*diag(lambda(:,t).^(1/4)); Cs = Cs*diag(lambda(:,t).^(1/4));Ds = Ds*diag(lambda(:,t).^(1/4));      
    
    Atemp = sparse(I,F); Btemp = sparse(J,F);Ctemp = sparse(K,F);Dtemp = sparse(L,F);
    Atemp(idx_i,:) = As; Btemp(idx_j,:) = Bs; Ctemp(idx_k,:) = Cs; Dtemp(idx_l,:) = Ds;
    
    lamb = zeros(F,1);
    %now, normalize the common part of every column
        for f = 1:F
           norm_a = norm(Atemp(fixed_i,f),2);Atemp(:,f) = Atemp(:,f)/norm_a; lamb(f) = norm_a;
           norm_b = norm(Btemp(fixed_j,f),2);Btemp(:,f) = Btemp(:,f)/norm_b; lamb(f) = lamb(f)*norm_b;
           norm_c = norm(Ctemp(fixed_k,f),2);Ctemp(:,f) = Ctemp(:,f)/norm_c; lamb(f) = lamb(f)*norm_c;
           norm_d = norm(Dtemp(fixed_l,f),2);Dtemp(:,f) = Dtemp(:,f)/norm_d; lamb(f) = lamb(f)*norm_d;
        end
        
%         save(sprintf('part%d',t),Atemp,Btemp,Ctemp,lamb);
        [ii jj kk] = find(Atemp);data_dump = [ii jj kk];dlmwrite(sprintf('temp/A%d',t),data_dump,'\t');
        [ii jj kk] = find(Btemp);data_dump = [ii jj kk];dlmwrite(sprintf('temp/B%d',t),data_dump,'\t');
        [ii jj kk] = find(Ctemp);data_dump = [ii jj kk];dlmwrite(sprintf('temp/C%d',t),data_dump,'\t');
        [ii jj kk] = find(Dtemp);data_dump = [ii jj kk];dlmwrite(sprintf('temp/D%d',t),data_dump,'\t');
        [ii jj kk] = find(lamb);data_dump = [ii jj kk];dlmwrite(sprintf('temp/lamb%d',t),data_dump,'\t');        

end

for t = 1:times
    
    data_dump = dlmread(sprintf('temp/A%d',t),'\t');Atemp2 = sparse(data_dump(:,1),data_dump(:,2),data_dump(:,3),size(A,1),size(A,2));
    data_dump = dlmread(sprintf('temp/B%d',t),'\t');Btemp2 = sparse(data_dump(:,1),data_dump(:,2),data_dump(:,3),size(B,1),size(B,2));
    data_dump = dlmread(sprintf('temp/C%d',t),'\t');Ctemp2 = sparse(data_dump(:,1),data_dump(:,2),data_dump(:,3),size(C,1),size(C,2));
    data_dump = dlmread(sprintf('temp/D%d',t),'\t');Dtemp2 = sparse(data_dump(:,1),data_dump(:,2),data_dump(:,3),size(D,1),size(D,2));
    data_dump = dlmread(sprintf('temp/lamb%d',t),'\t');lamb2 = spconvert(data_dump);





    lambda(:,t) = lamb2;
    %now, possibly reorder the updates, according to their similarity with
    %the existing vectors:
    if t ==1
        A = Atemp2; B = Btemp2; C = Ctemp2; D = Dtemp2;
        [ii jj kk] = find(Atemp2);
        dlmwrite(sprintf('temp/A-cAll-it%d',t),[ii jj kk],'delimiter','\t');
        [ii jj kk] = find(Btemp2);
        dlmwrite(sprintf('temp/B-cAll-it%d',t),[ii jj kk],'delimiter','\t');
        [ii jj kk] = find(Ctemp2);
        dlmwrite(sprintf('temp/C-cAll-it%d',t),[ii jj kk],'delimiter','\t');
        [ii jj kk] = find(Dtemp2);
        dlmwrite(sprintf('temp/D-cAll-it%d',t),[ii jj kk],'delimiter','\t');
    else
        valA=zeros(F,1);valB=zeros(F,1);valC=zeros(F,1);valD=zeros(F,1);
        allCorrespondencesA = zeros(F,1);allCorrespondencesB = zeros(F,1);allCorrespondencesC = zeros(F,1);allCorrespondencesD = zeros(F,1);

        usedA = [];usedB = []; usedC=[]; usedD = [];
        for f1 = 1:F
            for f2 = 1:F
                valA(f2) = Atemp2(fixed_i,f1)'* A(fixed_i,f2);%/(norm(Atemp(idx_i,k),2)*norm(A(idx_i,k),2));
                valB(f2) = Btemp2(fixed_j,f1)'* B(fixed_j,f2);%/(norm(Btemp(idx_j,k),2)*norm(B(idx_j,k),2));
                valC(f2) = Ctemp2(fixed_k,f1)'* C(fixed_k,f2);%/(norm(Ctemp(idx_k,k),2)*norm(C(idx_k,k),2));
                valD(f2) = Dtemp2(fixed_l,f1)'* D(fixed_l,f2);%/(norm(Ctemp(idx_k,k),2)*norm(C(idx_k,k),2));
            end
            [junk idx] = sort(valA,'descend');            
            ii = 1;while (find(idx(ii) == usedA))ii = ii +1; end;idx = idx(ii);usedA = [usedA idx];%use the highest non-already used matching column            allCorrespondencesA(f1) = idx;
% %             mask2 = find(A(:,idx)==0);
% % %             A(mask2,idx) = A(mask2,idx) + Atemp2(mask2,f1);%Update ONLY the zero values
% %             dlmwrite(sprintf('temp/A-c%d-it%d',f1,t),[mask2 f1*ones(size(mask2)) full(Atemp2(mask2,f1))],'delimiter','\t');
            allCorrespondencesA(f1) = idx;

            [junk idx] = sort(valB,'descend');            
            ii = 1;while (find(idx(ii) == usedB))ii = ii +1; end;idx = idx(ii);usedB = [usedB idx];
% %             mask2 = find(B(:,idx)==0);
% % %             B(mask2,idx) = B(mask2,idx) + Btemp2(mask2,f1);%Update ONLY the zero values
% %             dlmwrite(sprintf('temp/B-c%d-it%d',f1,t),[mask2 f1*ones(size(mask2)) full(Btemp2(mask2,f1))],'delimiter','\t');
            allCorrespondencesB(f1) = idx;
            
           [junk idx] = sort(valC,'descend');            
            ii = 1;while (find(idx(ii) == usedC))ii = ii +1; end;idx = idx(ii);usedC = [usedC idx];
% %             mask2 = find(C(:,idx)==0);
% % %             C(mask2,idx) = C(mask2,idx) + Ctemp2(mask2,f1);%Update ONLY the zero values
% %             dlmwrite(sprintf('temp/C-c%d-it%d',f1,t),[mask2 f1*ones(size(mask2)) full(Ctemp2(mask2,f1))],'delimiter','\t');
            allCorrespondencesC(f1) = idx;
            
           [junk idx] = sort(valD,'descend');            
            ii = 1;while (find(idx(ii) == usedD))ii = ii +1; end;idx = idx(ii);usedD = [usedD idx];
% %             mask2 = find(C(:,idx)==0);
% % %             C(mask2,idx) = C(mask2,idx) + Ctemp2(mask2,f1);%Update ONLY the zero values
% %             dlmwrite(sprintf('temp/C-c%d-it%d',f1,t),[mask2 f1*ones(size(mask2)) full(Ctemp2(mask2,f1))],'delimiter','\t');
            allCorrespondencesD(f1) = idx;

        end
        dlmwrite(sprintf('temp/corrA%d',t),allCorrespondencesA-1,'\t');%Java indexing
        dlmwrite(sprintf('temp/corrB%d',t),allCorrespondencesB-1,'\t');
        dlmwrite(sprintf('temp/corrC%d',t),allCorrespondencesC-1,'\t');
        dlmwrite(sprintf('temp/corrD%d',t),allCorrespondencesD-1,'\t');
    end
%     time_phase_two = time_phase_two + toc;
end

arg1 = sprintf('java  -Xmx%dm -jar factor-merge.jar temp/A temp/corrA %d %d',MAXMEM,F,times);
arg2 = sprintf('java  -Xmx%dm -jar factor-merge.jar temp/B temp/corrB %d %d',MAXMEM,F,times);
arg3 = sprintf('java  -Xmx%dm -jar factor-merge.jar temp/C temp/corrC %d %d',MAXMEM,F,times);
arg4 = sprintf('java  -Xmx%dm -jar factor-merge.jar temp/D temp/corrD %d %d',MAXMEM,F,times);

%choose whether we should do the factor-merge in parallel or just serially,
%based on the size of the tensor
if (prod(dims)>=10^15)%too large, just do it serially
    disp('Serial FactorMerge')
    system(arg1);
    system(arg2);
    system(arg3);
    system(arg4);
else %we can do it in parallel!
    disp('Parallel FactorMerge')
    cluster = parcluster;
    job = createJob(cluster);
    createTask(job,@system,1,{{arg1} {arg2} {arg3} {arg4}});
    submit(job);
    wait(job);
    out = fetchOutputs(job);
    delete(job);
end


A = load('temp/A-final');
A = spconvert(A);
if size(A,1) ~= I || size(A,2) ~= F
    A(I,F) = 0;
end

B = load('temp/B-final');
B = spconvert(B);
if size(B,1) ~= J || size(B,2) ~= F
    B(J,F) = 0;
end
C = load('temp/C-final');
C = spconvert(C);
if size(C,1) ~= K || size(C,2) ~= F
    C(K,F) = 0;
end

D = load('temp/D-final');
D = spconvert(D);
if size(D,1) ~= L || size(D,2) ~= F
    D(L,F) = 0;
end


% system('rm -r temp');
rmdir('temp','s');

system(sprintf('rm %s-*',FILENAME));
system(sprintf('rm %s1.txt*',FILENAME));
system(sprintf('rm %s2.txt*',FILENAME));
system(sprintf('rm %s3.txt*',FILENAME));
system(sprintf('rm %s4.txt*',FILENAME));