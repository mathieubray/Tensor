%Vagelis Papalexakis, 2012-2013
%School of Computer Science, Carnegie Mellon University
%Demonstrating the use of the 4-mode ParCube algorithm, using randomly generated
%tensors.
clear all;close all;clc
I = 100; J = 100; K = 100;L = 100;
F = 2;
density = 0.05;
Ao = sprand(I,F,density); Bo = sprand(J,F,density); Co = sprand(K,F,density); Do = sprand(L,F,density);

X = khatrirao(Co,khatrirao(Bo,Ao))*Do';
X = sptensor(reshape(X, [I J K L]));

ii = find(X);
FILENAME = 'test-data/test-4way';
dlmwrite(FILENAME,[ii X(ii)],'\t');

%this is how you call the function
matlabpool open %open a pool of parallel workers
sample_factor = [2 2 2 2];%sampling factor *per mode*
times = 5;%number of sampling repetitions. Ideally is at least equal to the number of cores (so that it runs entirely in parallel)
dims = [I J K L];
%F is already specified while creating the tensor. In general, F is the
%number of components to be extracted.
[A B C D lambda] = file_parCube_4mode(FILENAME,dims,F,sample_factor,times);
matlabpool close