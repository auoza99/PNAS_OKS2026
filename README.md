# PNAS_OKS2026
Codes for paper "Traveling waves in a continuum model of schooling swimmers" by A. U. Oza, E. Kanso and M. J. Shelley, published in the Proceedings of the National Academy of Sciences (PNAS 2026)

The simulation data in Figure 3 and Figure 6 can be obtained by running the MATLAB code ContinuumPDESolver.m, which solves the continuum PDE model (Eq. 6 in the paper). The initial data are provided as csv files in the appropriate folders. For Figure 3, the initial data are in the folder Fig3. The command for generating the data is

[tdata,rho_data,C_data,S_data,xp] = ContinuumPDESolver(readmatrix('rho_initial.csv'),readmatrix('C_initial.csv'),readmatrix('S_initial.csv'),3,0.05,5,512,2^(-8),5000);

For Figure 6, the initial data are in the folder Fig6. The command for generating the data is

[tdata,rho_data,C_data,S_data,xp] = ContinuumPDESolver(readmatrix('rho_initial.csv'),readmatrix('C_initial.csv'),readmatrix('S_initial.csv'),3,0.05,320,16384,2^(-12),800); toc;

For the quasistatic 



