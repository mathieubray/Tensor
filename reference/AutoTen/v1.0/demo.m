clear all;close all;clc
%Vagelis Papalexakis - Carnegie Mellon University, School of Computer
%Science (2014-2015)

allFo = 3:4;
maxiter = 2;
est_rank = zeros(length(allFo),maxiter);
est_rank_baseline1 = zeros(length(allFo),maxiter);
est_rank_baseline2 = zeros(length(allFo),maxiter);
I = 50;

for it = 1:maxiter
for Fo = allFo
    %for sparse data, use the following line
    info = create_problem('Factor_Generator',@randi,'Size',[I I I],'Num_Factors',Fo,'Lambda_Generator', @ones,'Sparse_Generation',500);%sparse
    %for dense data, use the following line
%     info = create_problem('Factor_Generator',@randi,'Size',[I I I],'Num_Factors',Fo,'Lambda_Generator', @ones);%dense

    Fmax = 2*Fo;
    [Fac, c, F_est,loss] = AutoTen(info.Data,Fmax,2);
    [Fac, F_est_baseline3] = AutoTenBaseline(info.Data,Fmax,1);
    [Fac, F_est_baseline3] = AutoTenBaseline(info.Data,Fmax,2);
    est_rank(find(Fo==allFo),it) = F_est;
    est_rank_baseline1(find(Fo==allFo),it) = F_est_baseline3;
    est_rank_baseline2(find(Fo==allFo),it) = F_est_baseline3;

end
end


a = abs(mean(repmat(allFo',1,maxiter) - est_rank,2));
b = abs(mean(repmat(allFo',1,maxiter) - est_rank_baseline1,2));
c = abs(mean(repmat(allFo',1,maxiter) - est_rank_baseline2,2));

figure
bar(allFo,[a b c])
xlabel('Rank');ylabel('Mean Error')
legend('AutoTen','Baseline_2','Baseline_3','Location','Best')