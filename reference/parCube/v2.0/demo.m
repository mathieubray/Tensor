%Vagelis Papalexakis, 2012
%School of Computer Science, Carnegie Mellon University
%Demonstrating the use of the ParCube algorithm, using randomly generated
%tensors.
clear all;close all;clc

%Memory based version
I = 100; J = 100; K = 100;
F = 5;
density = 0.05;
Ao = sprand(I,F,density); Bo = sprand(J,F,density); Co = sprand(K,F,density);

X = 10*khatrirao(Bo,Ao)*Co';
X = sptensor(reshape(X, [I J K]));

sample_factor = [2 2 2];
repetitions = 4;
nonneg = 1;
[A B C lambda] = parCube(X,F,sample_factor,repetitions,nonneg);


%Parallel, disk based version
times = 10;
FILENAME = 'test-data/test';
I = 101; J = 102; K = 103;
sample_factor = [2 2 2];
F = 4;
dims = [I J K];
matlabpool open %open a pool of parallel workers
[A B C lambda] = file_parCube(FILENAME,dims,F,sample_factor,times);
matlabpool close