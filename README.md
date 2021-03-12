# MTP2
## Simulation code for provenance recovery
### Project Overview : 
In WSNs, packets travel from source to destination node in a multi hop manner. 
Provenance records the data source, forwarding, and aggregating information of data packets on their way to the base station.
Provenance is critical for assessing the trustworthiness of the received data, diagnosing network failures, detecting early signs of attacks, etc.
![path](images/path.png)
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

