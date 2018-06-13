function [A B C lambda] = file_parCube(FILENAME,dims,F,sample_factor,times)
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
else
    s1 = sample_factor; 
    s2 = sample_factor;
    s3 = sample_factor;
end
I = dims(1); J = dims(2); K = dims(3);

MAXMEM = 8192;

Xma = sparse(I,1); Xmb = sparse(J,1); Xmc = sparse(K,1);
% tic
disp('calling external jar file')
system(sprintf('java  -Xmx%dm -jar tensor-sampling.jar %s %d %d %d',MAXMEM,FILENAME,I,J,K));
% [data] = dlmread(sprintf('%s1.txt',FILENAME),'\t');Xma(uint32(data(:,1))+1) = data(:,2);
% [data] = dlmread(sprintf('%s2.txt',FILENAME),'\t');Xmb(uint32(data(:,1))+1) = data(:,2);
% [data] = dlmread(sprintf('%s3.txt',FILENAME),'\t');Xmc(uint32(data(:,1))+1) = data(:,2);


FILENAME_IDX1 = sprintf('%s1.txt',FILENAME);
FILENAME_IDX2 = sprintf('%s2.txt',FILENAME);
FILENAME_IDX3 = sprintf('%s3.txt',FILENAME);

system('mkdir temp')

disp('Created the marginals!')



% time_phase_one = toc;
% time_phase_two = 0;

A = sparse(I,F); B = sparse(J,F);C = sparse(K,F);
lambda = zeros(F,times);

p = 0.35;
%initialize the fixed set of indices
% idx_i = find(Xma~=0);
% idx_i = randsample_weighted(idx_i,I/s1,full(Xma(idx_i)));
system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d 0',MAXMEM,FILENAME_IDX1,ceil(I/s1)));
idx_i = load(sprintf('%s.sample-0',FILENAME_IDX1));


% idx_j = find(Xmb~=0);
% idx_j = randsample_weighted(idx_j,J/s2,full(Xmb(idx_j)));
system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d 0 ',MAXMEM,FILENAME_IDX2,ceil(J/s2)));
idx_j = load(sprintf('%s.sample-0',FILENAME_IDX2));


% idx_k = find(Xmc~=0);
% idx_k = randsample_weighted(idx_k,K/s3,full(Xmc(idx_k)));
system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d 0',MAXMEM,FILENAME_IDX3,ceil(K/s3)));
idx_k = load(sprintf('%s.sample-0',FILENAME_IDX3));


fixed_i = idx_i(1:ceil(p*length(idx_i)));
fixed_j = idx_j(1:ceil(p*length(idx_j)));
fixed_k = idx_k(1:ceil(p*length(idx_k)));


parfor t = 1:times
disp(sprintf('Pass #%d',t))
    %now make one pass to the file and select only specific indices
    
% %     idx_i = find(Xma~=0);
% %     idx_i = randsample_weighted(idx_i,I/s1,full(Xma(idx_i)));
% % 
% %     idx_j = find(Xmb~=0);
% %     idx_j = randsample_weighted(idx_j,J/s2,full(Xmb(idx_j)));
% % 
% %     idx_k = find(Xmc~=0);
% %     idx_k = randsample_weighted(idx_k,K/s3,full(Xmc(idx_k)));

    system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d %d',MAXMEM,FILENAME_IDX1,ceil(I/s1),t));
    idx_i = load(sprintf('%s.sample-%d',FILENAME_IDX1,t));
    system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d %d',MAXMEM,FILENAME_IDX2,ceil(J/s2),t));
    idx_j = load(sprintf('%s.sample-%d',FILENAME_IDX2,t));
    system(sprintf('java  -Xmx%dm -jar randsample-weighted.jar %s %d %d',MAXMEM,FILENAME_IDX3,ceil(K/s3),t));
    idx_k = load(sprintf('%s.sample-%d',FILENAME_IDX3,t));
    
    
    idx_i = union(idx_i,fixed_i);dlmwrite(sprintf('%s.sample-%d',FILENAME_IDX1,t),idx_i,'delimiter','','precision','%.0f');
    idx_j = union(idx_j,fixed_j);dlmwrite(sprintf('%s.sample-%d',FILENAME_IDX2,t),idx_j,'delimiter','','precision','%.0f');
    idx_k = union(idx_k,fixed_k);dlmwrite(sprintf('%s.sample-%d',FILENAME_IDX3,t),idx_k,'delimiter','','precision','%.0f');
    
    system(sprintf('java  -Xmx%dm -jar create-sampled-tensor.jar %s.sample-%d %s.sample-%d %s.sample-%d %s %d',MAXMEM,FILENAME_IDX1,t,FILENAME_IDX2,t,FILENAME_IDX3,t,FILENAME,t));
    dat = dlmread(sprintf('%s-%d',FILENAME,t),'\t');
    
    [junk tmp_i] = ismember(dat(:,1),idx_i);
    [junk tmp_j] = ismember(dat(:,2),idx_j);
    [junk tmp_k] = ismember(dat(:,3),idx_k);
    
    Xs = sptensor([tmp_i tmp_j tmp_k],dat(:,end),[length(idx_i) length(idx_j) length(idx_k)]);   

    
% % %     Xs = sptensor([],[],[length(idx_i) length(idx_j) length(idx_k)]);
% % %     fid = fopen(FILENAME,'r');
% % %     current_count = 1;
% % %     line = fgets(fid);
% % %     while line~=-1
% % %         data = textscan(line,'%s','delimiter','\t');
% % %         data = data{1};
% % %         if length(data)~=4
% % %             warning(sprintf('Input row was ill formated: %s',line))
% % %             line = fgets(fid);
% % %             continue;
% % %         end
% % %         i = str2num(char(data(1)));
% % %         j = str2num(char(data(2)));
% % %         k = str2num(char(data(3)));
% % %         v = data(4); v = str2double(char(v));
% % %         cond1 = find(i == idx_i);
% % %         cond2 = find(j == idx_j);
% % %         cond3 = find(k == idx_k);
% % % 
% % %         if ~isempty(cond1) && ~isempty(cond2) && ~isempty(cond3)
% % %             Xs (cond1,cond2,cond3) = v;
% % %         end
% % %         line = fgets(fid);
% % %         
% % % 
% % %         current_count =current_count + 1;
% % %     end
% % %     fclose(fid);
    

%     tic
    factors = cp_nmu(Xs,F);
    As=sparse(factors.U{1});Bs=sparse(factors.U{2});Cs=sparse(factors.U{3});%lambda = factors.lambda;As = As*diag(lambda);
    lambda(:,t) = factors.lambda;
    As = As*diag(lambda(:,t).^(1/3)); Bs = Bs*diag(lambda(:,t).^(1/3)); Cs = Cs*diag(lambda(:,t).^(1/3));      
    
    Atemp = sparse(I,F); Btemp = sparse(J,F);Ctemp = sparse(K,F);
    Atemp(idx_i,:) = As; Btemp(idx_j,:) = Bs; Ctemp(idx_k,:) = Cs;
    
    lamb = zeros(F,1);
    %now, normalize the common part of every column
        for f = 1:F
           norm_a = norm(Atemp(fixed_i,f),2);Atemp(:,f) = Atemp(:,f)/(norm_a+eps); lamb(f) = norm_a;
           norm_b = norm(Btemp(fixed_j,f),2);Btemp(:,f) = Btemp(:,f)/(norm_b+eps); lamb(f) = lamb(f)*norm_b;
           norm_c = norm(Ctemp(fixed_k,f),2);Ctemp(:,f) = Ctemp(:,f)/(norm_c+eps); lamb(f) = lamb(f)*norm_c;
        end
        
%         save(sprintf('part%d',t),Atemp,Btemp,Ctemp,lamb);
        [ii jj kk] = find(Atemp);data_dump = [ii jj kk];dlmwrite(sprintf('temp/A%d',t),data_dump,'\t');
        [ii jj kk] = find(Btemp);data_dump = [ii jj kk];dlmwrite(sprintf('temp/B%d',t),data_dump,'\t');
        [ii jj kk] = find(Ctemp);data_dump = [ii jj kk];dlmwrite(sprintf('temp/C%d',t),data_dump,'\t');
        [ii jj kk] = find(lamb);data_dump = [ii jj kk];dlmwrite(sprintf('temp/lamb%d',t),data_dump,'\t');
        

end

for t = 1:times
    
    data_dump = dlmread(sprintf('temp/A%d',t),'\t');Atemp2 = sparse(data_dump(:,1),data_dump(:,2),data_dump(:,3),size(A,1),size(A,2));
    data_dump = dlmread(sprintf('temp/B%d',t),'\t');Btemp2 = sparse(data_dump(:,1),data_dump(:,2),data_dump(:,3),size(B,1),size(B,2));
    data_dump = dlmread(sprintf('temp/C%d',t),'\t');Ctemp2 = sparse(data_dump(:,1),data_dump(:,2),data_dump(:,3),size(C,1),size(C,2));
    data_dump = dlmread(sprintf('temp/lamb%d',t),'\t');lamb2 = spconvert(data_dump);





    lambda(:,t) = lamb2;
    %now, possibly reorder the updates, according to their similarity with
    %the existing vectors:
    if t ==1
        t
        A = Atemp2; B = Btemp2; C = Ctemp2;
        [ii jj kk] = find(Atemp2);
        dlmwrite(sprintf('temp/A-cAll-it%d',t),[ii jj kk],'delimiter','\t');
        [ii jj kk] = find(Btemp2);
        dlmwrite(sprintf('temp/B-cAll-it%d',t),[ii jj kk],'delimiter','\t');
        [ii jj kk] = find(Ctemp2);
        dlmwrite(sprintf('temp/C-cAll-it%d',t),[ii jj kk],'delimiter','\t');
    else
        t
        valA=zeros(F,1);valB=zeros(F,1);valC=zeros(F,1);
        allCorrespondencesA = zeros(F,1);allCorrespondencesB = zeros(F,1);allCorrespondencesC = zeros(F,1);

        usedA = [];usedB = []; usedC=[];
        for f1 = 1:F
            for f2 = 1:F
                valA(f2) = Atemp2(fixed_i,f1)'* A(fixed_i,f2);%/(norm(Atemp(idx_i,k),2)*norm(A(idx_i,k),2));
                valB(f2) = Btemp2(fixed_j,f1)'* B(fixed_j,f2);%/(norm(Btemp(idx_j,k),2)*norm(B(idx_j,k),2));
                valC(f2) = Ctemp2(fixed_k,f1)'* C(fixed_k,f2);%/(norm(Ctemp(idx_k,k),2)*norm(C(idx_k,k),2));
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
        end
        dlmwrite(sprintf('temp/corrA%d',t),allCorrespondencesA-1,'\t');%Java indexing
        dlmwrite(sprintf('temp/corrB%d',t),allCorrespondencesB-1,'\t');
        dlmwrite(sprintf('temp/corrC%d',t),allCorrespondencesC-1,'\t');
    end
%     time_phase_two = time_phase_two + toc;
% disp(sprintf('Merged %d out of %d parts',t,times))

end
arg1 = sprintf('java  -Xmx%dm -jar factor-merge.jar temp/A temp/corrA %d %d',MAXMEM,F,times);
arg2 = sprintf('java  -Xmx%dm -jar factor-merge.jar temp/B temp/corrB %d %d',MAXMEM,F,times);
arg3 = sprintf('java  -Xmx%dm -jar factor-merge.jar temp/C temp/corrC %d %d',MAXMEM,F,times);

%choose whether we should do the factor-merge in parallel or just serially,
%based on the size of the tensor
if (prod(dims)>=10^15)%too large, just do it serially
    disp('Serial FactorMerge')
    system(arg1);
    system(arg2);
    system(arg3);
else %we can do it in parallel!
    disp('Parallel FactorMerge')
    cluster = parcluster;
    job = createJob(cluster);
    createTask(job,@system,1,{{arg1} {arg2} {arg3}});
    submit(job);
    wait(job);
    out = fetchOutputs(job);
    delete(job);
end
% system(sprintf('java  -Xmx%dm -jar factor-merge.jar temp/A temp/corrA %d %d &',MAXMEM,F,times));
% disp('Done with A')
% system(sprintf('java  -Xmx%dm -jar factor-merge.jar temp/B temp/corrB %d %d ',MAXMEM,F,times));
% disp('Done with B')
% system(sprintf('java  -Xmx%dm -jar factor-merge.jar temp/C temp/corrC %d %d ',MAXMEM,F,times));
% disp('Done with C')
% 
% %NEED SYNCHRONIZATION POINT HERE AND READ DATA TO VARS A B C
% count=0;
% while(~exist('temp/A-final') && ~exist('temp/B-final'))%block until files are ready
%     count = count+1;
% end
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

rmdir('temp','s');

% total_time = time_phase_one + time_phase_two;
% disp(sprintf('Time elapsed: %d seconds',total_time))

lambda = mean(lambda,2);
% lambda = max(lambda');
% A = A*diag(lambda);
% % A(A<10^-8) = 0;
% % B(B<10^-8) = 0;
% % C(C<10^-8) = 0;

system(sprintf('rm %s-*',FILENAME));
system(sprintf('rm %s1.txt*',FILENAME));
system(sprintf('rm %s2.txt*',FILENAME));
system(sprintf('rm %s3.txt*',FILENAME));