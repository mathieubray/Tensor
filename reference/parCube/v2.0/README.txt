File demo.m contains an example of executing
	-parCube.m: memory based implementation of ParCube
	-file_parCube.m: parallel, disk based implementation of ParCube (this file, currently runs on UNIX based systems)
File demo_4mode.m contains an example of executing the file based ParCube for 4-mode tensors.

[June 2013 Update]
Most of the file based functions are now implemented in Java and the code uses the compiled .jar files.

If you use this code, please cite our paper:

@incollection{
year={2012},
isbn={978-3-642-33459-7},
booktitle={Machine Learning and Knowledge Discovery in Databases},
volume={7523},
series={Lecture Notes in Computer Science},
editor={Flach, PeterA. and Bie, Tijl and Cristianini, Nello},
doi={10.1007/978-3-642-33460-3_39},
title={ParCube: Sparse Parallelizable Tensor Decompositions},
url={http://dx.doi.org/10.1007/978-3-642-33460-3_39},
publisher={Springer Berlin Heidelberg},
keywords={Tensors; PARAFAC decomposition; sparsity; sampling; randomized algorithms; parallel algorithms},
author={Papalexakis, EvangelosE. and Faloutsos, Christos and Sidiropoulos, NicholasD.},
pages={521-536}
}
