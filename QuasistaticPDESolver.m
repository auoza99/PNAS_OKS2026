function [tdata,rho_data] = QuasistaticPDESolver(rho0,alpha,nu,L,N,dt,tmax)

% This code solves the quasistatic PDE model (Eq. 14 in paper) using the 
% ... numerical method described in Section 4 of the SI Appendix. We use a
% ... pseudospectral method in space and a fourth-order Runge-Kutta
% ... time-stepping scheme with an integrating factor.

% INPUTS:
% rho0: initial density rho (rho(x,0)), a vector of length N
% alpha: model parameter (dimensionless flapping frequency)
% nu: model parameter (dimensionless diffusivity)
% L: sets size of spatial domain, which is [-L*pi,L*pi]
% N: number of grid points in space
% dt: time step
% tmax: final time to run the simulation up to

% OUTPUTS:
% tdata: vector of time values
% rho_data: matrix consisting of the solution rho(x,t): row i corresponds
% ... to rho(x,tdata(i))

% numerical simulation parameters
p = -round(log2(dt)); % time step is dt = 2^(-p)
k = [0:N/2-1 0 -N/2+1:-1]; % k-vector
k2 = [0:N/2-1 1*N/2 -N/2+1:-1]; 
x = 2*pi*L/N*(-N/2:N/2-1); % x-vector
Mval = nu*k2.^2/L^2 - 1i*k/L;

rho_hat = fft(rho0); % Fourier transform of rho(x,0)

EM = exp(-Mval*dt);
EM2 = exp(-Mval*dt/2);

nplt = ceil(1/dt); % sets number of time steps to store
nmax = round(tmax/dt); % total number of time steps
rho_data = NaN*ones(ceil(nmax/nplt)+1,N); % store \rho
rho_data(1,:) = rho0;

tdata = NaN*ones(1,ceil(nmax/nplt)+1); % store time values
tdata(1) = 0;

cnt = 2; % counter to store rho data

% time-stepping loop, given in Eqs. 48-49 of SI Appendix
for n = 1:nmax
    
    t = n*dt;
    
    % compute f1
    f1 = NFunc(rho_hat,alpha,nu,N,L);
    f1 = fft(f1); 
    
    % compute f2
    in_1 = EM2.*(rho_hat + dt/2*f1);
    f2 = NFunc(in_1,alpha,nu,N,L);
    f2 = fft(f2); 
    
    % compute f3
    in_1 = EM2.*rho_hat + dt/2*f2;
    f3 = NFunc(in_1,alpha,nu,N,L);
    f3 = fft(f3); 
    
    % compute f4
    in_1 = EM.*rho_hat + dt*EM2.*f3;
    f4 = NFunc(in_1,alpha,nu,N,L);
    f4 = fft(f4); 
    
    % perform time step
    rho_hat = EM.*rho_hat + dt/6*(EM.*f1 + 2*EM2.*(f2 + f3) + f4);
        
    % store \rho(x,t)
    if (mod(n,nplt) == 0)
        tdata(cnt) = t;
        
        rho = real(ifft(rho_hat));
        rho_data(cnt,:) = rho; 

        cnt = cnt + 1;
    end;
end;

function [out] = NFunc(rho_hat,alpha,nu,N,L) 
% nonlinear term in Eq. 14 written in real space

k = [0:N/2-1 0 -N/2+1:-1];

rho = real(ifft(rho_hat));
rhox = 1/L*real(ifft(1i*k.*rho_hat)); % partial derivative \partial_x\rho

out = -1/L*real(ifft(1i*k.*fft(rho./(1 + alpha^2 + rho).*(rho + nu*rhox))));