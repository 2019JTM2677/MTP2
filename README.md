# MTP2
## Simulation code for provenance recovery
### Project Overview : 
In WSNs, packets travel from source to destination node in a multi hop manner. 
Provenance records the data source, forwarding, and aggregating information of data packets on their way to the base station.
Provenance is critical for assessing the trustworthiness of the received data, diagnosing network failures, detecting early signs of attacks, etc.

![path](https://github.com/2019JTM2677/MTP2/blob/master/images/path.PNG)

*Path in n-node network with relaying mobile nodes* 

We use compressive sensing technique to embed the information about the realying nodes into the packet using edge embedding and duble edge embedding methods. 
CVX and OMP and its variants are used to solve the compressive sensing problem and the error rate in path recovery for these methods are compared.

### Prerequisites :
- MATLAB (R2020b)
- CVX optimisation toolbox
To install cvx toolbox in matlab, follow the steps:
Download the cvx from this link http://cvxr.com/cvx/download/. This page also mentions the installation method.
Download the file for your platform and then place it at a convenient location. Then from Matlab change the current folder to the cvx folder and run:
```
cvx_setup
```
Please note that the use MOSEK require a CVX Professional license key, for full information go through this page:
http://cvxr.com/cvx/doc/mosek.html

### How to run the code:
1. Install all the prerequisites that are mentioned above which is required to run the code
2. Set parameters n= no. of nodes in the network, h= hop length of path required, m= no. of measurements (column size of sensing matrix), error threshold, no of pkts for simulation

### Some results:
![SE DE CVX OMP](https://github.com/2019JTM2677/MTP2/blob/master/images/SE%20vs%20DE%20(cvx%2Comp).png)

*SE vs DE for CVX and OMP*

![SE OMP vs modified OMP](https://github.com/2019JTM2677/MTP2/blob/master/images/SE_OMP%20vs%20mod%20OMP.png)

*SE OMP vs modified OMP*

![DE OMP vs modified OMP](https://github.com/2019JTM2677/MTP2/blob/master/images/DE_OMP_vs_modOMP.png)

*DE OMP vs modified OMP*

![SE topo 1](https://github.com/2019JTM2677/MTP2/blob/master/images/topo%20n%3D7%20h%3D4.PNG)
![SE topo 2](https://github.com/2019JTM2677/MTP2/blob/master/images/topo5.PNG)

*SE OMP vs modified OMP vs modified OMP with topology knowledge*

![DE topo 1](https://github.com/2019JTM2677/MTP2/blob/master/images/de_topo1.PNG)
![DE topo 2](https://github.com/2019JTM2677/MTP2/blob/master/images/de_topo2.PNG)

*DE OMP vs modified OMP vs modified OMP with topology knowledge*
