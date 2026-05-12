function [tdata,rho_data,C_data,S_data,xp] = ContinuumPDESolver(rho0,C0,S0,alpha,nu,L,N,dt,tmax)

% This code solves the continuum PDE model (Eq. 6 in paper) using the 
% ... numerical method described in Section 4 of the SI Appendix. We use a
% ... pseudospectral method in space and a fourth-order Runge-Kutta
% ... time-stepping scheme with an integrating factor.

% INPUTS:
% rho0: initial density rho (rho(x,0)), a vector of length N
% C0: initial field C (C(x,0)), a vector of length N
% S0: initial field S (S(x,0)), a vector of length N
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
% C_data: matrix consisting of the solution C(x,t)
% S_data: matrix consisting of the solution S(x,t)
% xp: matrix consisting of "Lagrangian" particles displayed in Fig. 3
% ... of the paper. Each row has the trajectory of a different particle.

% numerical simulation parameters
p = -round(log2(dt)); % time step is dt = 2^(-p)
k = [0:N/2-1 0 -N/2+1:-1]; % k-vector
k2 = [0:N/2-1 1*N/2 -N/2+1:-1];
x = 2*pi*L/N*(-N/2:N/2-1); % x-vector
ks = k/L + 1i*nu*k2.^2/L^2; % derivatives in Fourier space
fc = (1 + 1i*nu*k/L)./((1 + 1i*ks + 1i*alpha).*(1 + 1i*ks - 1i*alpha)); % phi defined in Eq. 51 of SI Appendix

rho_hat = fft(rho0); % Fourier transform of rho(x,0)
C_hat = fft(C0); % Fourier transform of C(x,0)
S_hat = fft(S0); % Fourier transform of S(x,0)

% matrix elements of exp(-M*dt), given in Eq. 50 of SI Appendix
Mdt_11 = exp(1i*ks*dt);
Mdt_12 = zeros(1,N);
Mdt_13 = zeros(1,N);
Mdt_21 = fc.*((1 + 1i*ks).*(exp(-dt)*cos(alpha*dt)*ones(1,N) - exp(1i*ks*dt)) - alpha*exp(-dt)*sin(alpha*dt)*ones(1,N));
Mdt_22 = exp(-dt)*cos(alpha*dt)*ones(1,N);
Mdt_23 = -exp(-dt)*sin(alpha*dt)*ones(1,N);
Mdt_31 = fc.*(alpha*exp(-dt)*cos(alpha*dt)*ones(1,N) - alpha*exp(1i*ks*dt) + (1 + 1i*ks)*exp(-dt)*sin(alpha*dt));
Mdt_32 = exp(-dt)*sin(alpha*dt)*ones(1,N);
Mdt_33 = exp(-dt)*cos(alpha*dt)*ones(1,N);

% matrix elements of exp(-M*dt/2), given in Eq. 50 of SI Appendix
Mdt2_11 = exp(1i*ks*dt/2);
Mdt2_12 = zeros(1,N);
Mdt2_13 = zeros(1,N);
Mdt2_21 = fc.*((1 + 1i*ks).*(exp(-dt/2)*cos(alpha*dt/2)*ones(1,N) - exp(1i*ks*dt/2)) - alpha*exp(-dt/2)*sin(alpha*dt/2)*ones(1,N));
Mdt2_22 = exp(-dt/2)*cos(alpha*dt/2)*ones(1,N);
Mdt2_23 = -exp(-dt/2)*sin(alpha*dt/2)*ones(1,N);
Mdt2_31 = fc.*(alpha*exp(-dt/2)*cos(alpha*dt/2)*ones(1,N) - alpha*exp(1i*ks*dt/2) + (1 + 1i*ks)*exp(-dt/2)*sin(alpha*dt/2));
Mdt2_32 = exp(-dt/2)*sin(alpha*dt/2)*ones(1,N);
Mdt2_33 = exp(-dt/2)*cos(alpha*dt/2)*ones(1,N);

nplt = ceil(1/dt); % sets number of time steps to store
nmax = round(tmax/dt); % total number of time steps
rho_data = NaN*ones(ceil(nmax/nplt)+1,N); % store \rho
C_data = NaN*ones(ceil(nmax/nplt)+1,N); % store C
S_data = NaN*ones(ceil(nmax/nplt)+1,N); % store S
rho_data(1,:) = rho0;
C_data(1,:) = C0;
S_data(1,:) = S0;

tdata = NaN*ones(1,ceil(nmax/nplt)+1); % store time values
tdata(1) = 0;

% Generate Lagrangian particle initial positions with uniform distribution
nump = 50; % number of Lagrangian particles
xp = NaN*ones(nump,nmax); % position of Lagrangian particle
xuni = linspace(-L*pi,L*pi,nump+2);
xp(:,1) = xuni(2:end-1)';

cnt = 2; % counter to store rho, C and S data

% time-stepping loop, given in Eqs. 48-49 of SI Appendix
for n = 1:nmax
    
    t = n*dt; % time value
    
    % compute f1
    [f1_1,f1_2,f1_3] = NFunc(rho_hat,C_hat,N,L);
    f1_1 = fft(f1_1); 
    f1_2 = fft(f1_2); 
    f1_3 = fft(f1_3);
    
    % compute f2
    sum_1 = rho_hat + dt/2*f1_1;
    sum_2 = C_hat + dt/2*f1_2;
    sum_3 = S_hat + dt/2*f1_3;
    in_1 = Mdt2_11.*sum_1 + Mdt2_12.*sum_2 + Mdt2_13.*sum_3;
    in_2 = Mdt2_21.*sum_1 + Mdt2_22.*sum_2 + Mdt2_23.*sum_3;
    [f2_1,f2_2,f2_3] = NFunc(in_1,in_2,N,L);
    f2_1 = fft(f2_1); 
    f2_2 = fft(f2_2); 
    f2_3 = fft(f2_3);
    
    % compute f3
    in_1 = Mdt2_11.*rho_hat + Mdt2_12.*C_hat + Mdt2_13.*S_hat + dt/2*f2_1;
    in_2 = Mdt2_21.*rho_hat + Mdt2_22.*C_hat + Mdt2_23.*S_hat + dt/2*f2_2;
    [f3_1,f3_2,f3_3] = NFunc(in_1,in_2,N,L);
    f3_1 = fft(f3_1); 
    f3_2 = fft(f3_2); 
    f3_3 = fft(f3_3);
    
    % compute f4
    in_1a = Mdt_11.*rho_hat + Mdt_12.*C_hat + Mdt_13.*S_hat;
    in_2a = Mdt_21.*rho_hat + Mdt_22.*C_hat + Mdt_23.*S_hat;
    in_1b = Mdt2_11.*f3_1 + Mdt2_12.*f3_2 + Mdt2_13.*f3_3;
    in_2b = Mdt2_21.*f3_1 + Mdt2_22.*f3_2 + Mdt2_23.*f3_3;
    
    in_1 = in_1a + dt*in_1b;
    in_2 = in_2a + dt*in_2b;
    [f4_1,f4_2,f4_3] = NFunc(in_1,in_2,N,L);
    f4_1 = fft(f4_1); 
    f4_2 = fft(f4_2); 
    f4_3 = fft(f4_3);
    
    % compute zeroth term in time step
    s0_1 = Mdt_11.*rho_hat + Mdt_12.*C_hat + Mdt_13.*S_hat;
    s0_2 = Mdt_21.*rho_hat + Mdt_22.*C_hat + Mdt_23.*S_hat;
    s0_3 = Mdt_31.*rho_hat + Mdt_32.*C_hat + Mdt_33.*S_hat;
    
    % compute first term in time step
    s1_1 = Mdt_11.*f1_1 + Mdt_12.*f1_2 + Mdt_13.*f1_3;
    s1_2 = Mdt_21.*f1_1 + Mdt_22.*f1_2 + Mdt_23.*f1_3;
    s1_3 = Mdt_31.*f1_1 + Mdt_32.*f1_2 + Mdt_33.*f1_3;
    
    % compute second term in time step
    s2_1 = Mdt2_11.*f2_1 + Mdt2_12.*f2_2 + Mdt2_13.*f2_3;
    s2_2 = Mdt2_21.*f2_1 + Mdt2_22.*f2_2 + Mdt2_23.*f2_3;
    s2_3 = Mdt2_31.*f2_1 + Mdt2_32.*f2_2 + Mdt2_33.*f2_3;
    
    % compute third term in time step
    s3_1 = Mdt2_11.*f3_1 + Mdt2_12.*f3_2 + Mdt2_13.*f3_3;
    s3_2 = Mdt2_21.*f3_1 + Mdt2_22.*f3_2 + Mdt2_23.*f3_3;
    s3_3 = Mdt2_31.*f3_1 + Mdt2_32.*f3_2 + Mdt2_33.*f3_3;
    
    % perform time step
    rho_hat = s0_1 + dt/6*(s1_1 + 2*(s2_1 + s3_1) + f4_1);
    C_hat = s0_2 + dt/6*(s1_2 + 2*(s2_2 + s3_2) + f4_2);
    S_hat = s0_3 + dt/6*(s1_3 + 2*(s2_3 + s3_3) + f4_3);
    
    % evolve Lagrangian particle positions using forward Euler
    uvel = -1 - real(ifft(C_hat));
    for pind = 1:nump
        xp(pind,n+1) = xp(pind,n) + dt*interp1([x L*pi],[uvel uvel(1)],xp(pind,n));
        if (xp(pind,n+1) < -pi*L) % reflect particles periodically if they leave domain
            xp(pind,n+1) = xp(pind,n+1) + 2*pi*L;
        end;
    end;
        
    % store solutions \rho,C,S
    if (mod(n,nplt) == 0)
        tdata(cnt) = t;
        
        rho = real(ifft(rho_hat));
        C = real(ifft(C_hat));
        S = real(ifft(S_hat));

        rho_data(cnt,:) = rho;
        C_data(cnt,:) = C;
        S_data(cnt,:) = S;    

        cnt = cnt + 1;
    end;

end;

function [out1,out2,out3] = NFunc(rho_hat,C_hat,N,L)
% nonlinear terms in Eq. 6, written in real space

k = [0:N/2-1 0 -N/2+1:-1];

rho = real(ifft(rho_hat));
C = real(ifft(C_hat));

out1 = 1/L*real(ifft(1i*k.*fft(rho.*C))); % Eq. 6a
out2 = -rho.*C; % Eq. 6b
out3 = zeros(1,N); % Eq. 6c