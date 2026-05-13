function TravelingWaveBranches(alpha,nu,L,N,n)

% This code constructs branches of traveling wave solutions to the 
% ... continuum PDE model. It solves Eq. 12 in the paper using numerical
% ... continuation. The algorithm is described in Section 5 of the SI
% ... Appendix. 

% INPUTS:
% alpha: model parameter (dimensionless flapping frequency)
% nu: model parameter (dimensionless diffusivity)
% L: sets size of spatial domain, which is [-L*pi,L*pi]
% N: number of grid points in space
% n: index of the traveling wave branch (number of minima in the wave). The
% ... current settings produce all of the branches in Figures 4 and 5 of
% ... the paper except n = 1,5-7, for which different settings are needed.

% OUTPUTS:
% The density rho, fields C and S, and wave speed c are saved as variables
% ... R_out, C_out, S_out and c_out, respectively. The "exitflag", an
% ... output of MATLAB's root-finding algorithm, is also saved as "eflag".
% ... (Positive values of eflag indicate that the solver converged).
% ... The outputs are saved in .mat files, stored in folders named 
% ... n1, n2, etc. depending on the input n. The .mat files are labeled
% ... either A###.mat or R###.mat; the former (latter) is used if the 
% ... solution was obtained using continuation in the amplitude (density),
% ... with the value of the amplitude (density) in the file name.

options = optimset('display','off');

% First find values of wave speed c that demarcate the ends of the solution
% ... branch, as shown in SI Appendix Fig. S4d.
acc_c = 5; % number of digits of accuracy in c
dc = 10^(-acc_c);

% Solve Eq. 64 in SI Appendix
c1 = fsolve(@(c) c.^(3/2).*(1-sqrt(c))-nu*(alpha^2+1)/(alpha^2-1),.5); 
c2 = fsolve(@(c) c.^(3/2).*(1-sqrt(c))-nu*(alpha^2+1)/(alpha^2-1),.9);
cvec = linspace(c2,c1,(c2-c1)/dc+1); % range of possible c values
d_guess = -0.07; % initial guess for parameter d, in Eq. 53 of SI Appendix

dvec = NaN*ones(1,length(cvec)); % vector of d-values for each c-value

% eigenvalues of matrix M-tilde, defined in Eq. 57 of SI Appendix, as a
% ... function of c along the solution branch.
eval_vec = NaN*ones(3,length(cvec)); 

for ind = 1:length(cvec) % loop over c values
    c = cvec(ind);

    % find the value of d for which the largest eigenvalue of M-tilde
    % ... is purely imaginary.
    [d,~,eflag] = fsolve(@(x) func_dCalc(x,c,alpha,nu),d_guess,options);
    
    % store the value of d and the eigenvalue of M-tilde
    if (eflag > 0)
        d_guess = d; % update the guess if solver converged
        dvec(ind) = d;
        eval_vec(:,ind) = eig(Matrix(d,c,alpha,nu));
    end;
end;

% number of solution branches permitted over the periodic domain
nmax = floor(max(imag(eval_vec(2,:))*L)); 
if (n > nmax)
    disp('ERROR: n is too big')
end;

% store the endpoints of the curve in the (c,d) plane from which the
% ... solution branch bifurcates; see Fig. S4d in SI Appendix.
cEnds = NaN*[1 1];
dEnds = NaN*[1 1];

% for a given n, now find values of c and d for which the largest
% ... eigenvalue of M-tilde is 1i*n/L. Use the prior results to obtain a
% ... good initial guess.
ds = imag(eval_vec(2,:)) - n/L; 
loc = find(ds(2:end).*ds(1:end-1) < 0); % initial guesses for d
for lind = 1:2 
    [x_out,~,eflag] = fsolve(@(x) func(x,n/L,alpha,nu),[cvec(loc(lind)); dvec(loc(lind))],options); 
    cEnds(lind) = x_out(1);
    dEnds(lind) = x_out(2);
end;

% Now do numerical continuation, starting from the large-rho side of the
% ... solution branch. 
% Set parameters for MATLAB fsolve algorithm.
tol_fun = 1e-6; 
opt_tol = 1e-4;
options = optimoptions('fsolve','display','iter-detailed','MaxFunEvals',5000*(3*N+1),'FunctionTolerance',tol_fun,'OptimalityTolerance',opt_tol,'MaxIter',4000);

% Set initial guess for the solver
c = cEnds(2); 
d = dEnds(2); 
A_min = 0.01; % small amplitude perturbation from uniform base state
rts = roots([c -(d + (1 + alpha^2)*(1 - c)) -d*(1 + alpha^2)]); % Solve Eq. 55 in SI Appendix
R0 = rts(2); % mean density of the solution

% Initial guess according to Eq. 66 in SI Appendix
x = 2*pi*L/N*(-N/2:N/2-1);
Rg = R0 + A_min*sin(n*x/L); 
Cg = -Rg./(1 + alpha^2 + Rg);
Sg = alpha*Cg;

% Solve ODEs in Eq. 12 using a spectral method
[zg,~,eflag] = fsolve(@(z) TWaveFunc(z,alpha,nu,L,N,A_min,2),[Rg Cg Sg c],options); 

Rg = zg(1:N); % initial guess for density rho
Cg = zg(N+1:2*N); % initial guess for field C
Sg = zg(2*N+1:3*N); % initial guess for field S
cg = zg(3*N+1); % initial guess for wave speed c

% Perform numerical continuation in the wave amplitude A
Aval = 0.02; % starting amplitude
dA = 0.01; % increment in amplitude
eflag = Inf;
rho_prev = 1/N*trapz([Rg Rg(1)]);
drho = Inf;

% keep incrementing until solver fails (eflag criterion), or solution 
% ... branch becomes too steep (drho criterion)
while (eflag > 0 && ~(drho < -0.01)) 
    [z_out,~,eflag] = fsolve(@(z) TWaveFunc(z,alpha,nu,L,N,Aval,2),[Rg Cg Sg cg],options); 
    R_out = z_out(1:N);
    C_out = z_out(N+1:2*N);
    S_out = z_out(2*N+1:3*N);
    c_out = z_out(3*N+1);
    
    if (eflag > 0) % save outputs of solver if convergence is achieved
        save(['n' num2str(n) '/A' num2str(Aval) '.mat'],'R_out','C_out','S_out','c_out','eflag');
        Rg = R_out;
        Cg = C_out;
        Sg = S_out;
        cg = c_out;
        
        rho_current = 1/N*trapz([R_out R_out(1)]); % mean density of traveling wave solution

        disp(['A = ' num2str(Aval) '; rho = ' num2str(rho_current)]);
        
        drho = rho_current - rho_prev; % measure how much by which the mean density has changed
        rho_prev = rho_current;
        
        Aval = Aval + dA;
    else
        disp(['Break at A = ' num2str(Aval) '; eflag = ' num2str(eflag)])
        
        Aval = Aval + dA;
        break;
    end;
end;

% Now perform numerical continuation in the mean density
c = cEnds(1); % values of c and d at the small-rho end of branch
d = dEnds(1);
rts = roots([c -(d + (1 + alpha^2)*(1 - c)) -d*(1 + alpha^2)]); % Solve Eq. 55 in SI Appendix
rho_left = rts(2);
rho_val = rho_current;

% construct vector of mean densities to sweep over (rho_v)
acc_rho = 2; % number of digits of resolution in the vector of mean densities to sweep over
drho = 10^(-acc_rho);
rho_val_2 = round(rho_val,acc_rho);
if (rho_val_2 > rho_val)
    rho_max = rho_val_2 - drho;
else
    rho_max = rho_val_2;
end;

rho_min_2 = round(rho_left,acc_rho);
if (rho_min_2 > rho_left)
    rho_min = rho_min_2 - drho;
else
    rho_min = rho_min_2;
end;
rho_v = rho_max:-drho:rho_min;

% Loop over values of the mean density
for ind = 1:length(rho_v)
    m_prev = 1/N*trapz([Rg Rg(1)]); % value of the mean density
    [z_out,~,eflag] = fsolve(@(z) TWaveFunc(z,alpha,nu,L,N,rho_v(ind),1),[rho_v(ind)/m_prev*Rg Cg Sg cg],options); 
    R_out = z_out(1:N);
    C_out = z_out(N+1:2*N);
    S_out = z_out(2*N+1:3*N);
    c_out = z_out(3*N+1);

    if (eflag > 0) % save outputs of solver if convergence is achieved
        save(['n' num2str(n) '/R' num2str(rho_v(ind)) '.mat'],'R_out','C_out','S_out','c_out','eflag');
        Rg = R_out;
        Cg = C_out;
        Sg = S_out;
        cg = c_out;

        disp(['rho = ' num2str(rho_v(ind)) '; A = ' num2str(sqrt(2/N*trapz(([R_out R_out(1)]-rho_v(ind)).^2)))]);
    else
        disp(['Break at rho = ' num2str(rho_v(ind)) '; eflag = ' num2str(eflag)])
        break;
    end;
end;


function out = TWaveFunc(zv,alpha,nu,L,N,param,tpe)
% function that defines the ODEs in Eq. 12, satisfied by traveling wave
% ... solutions.  
 
% INPUTS:
% zv:  vector of length 3*N+1 that contains the unknowns: N each for rho,
% ... C and S, respectively; and one for c.
% alpha, nu, L and N are parameters defined at beginning of the code.
% param: mean value of density sought (R_int) for tpe = 1; 
% ... wave amplitude sought (Amp) for tpe = 2.
% tpe: determines whether we are finding solutions with a particular 
% ... mean density (tpe = 1) or amplitude (tpe = 2).

% OUTPUT: residual of the ODEs in Eq. 12.

R = zv(1:N); % density rho
C = zv(N+1:2*N); % field C
S = zv(2*N+1:3*N); % field S
c = zv(3*N+1); % wave speed c

k = [0:N/2-1 0 -N/2+1:-1]/L; % wavenumber k
k2 = [0:N/2-1 N/2 -N/2+1:-1]/L;

% Eq. 12 is computed using a spectral method
eq1 = -nu*real(ifft(-(k2.^2).*fft(R))) + c*real(ifft(1i*k.*fft(R))) + real(ifft(1i*k.*fft((-ones(1,N) - C).*R)));
eq2 = c*real(ifft(1i*k.*fft(C))) - R.*(-ones(1,N) - C) + nu*real(ifft(1i*k.*fft(R))) + C + alpha*S;
eq3 = c*real(ifft(1i*k.*fft(S))) + S - alpha*C;

R_int = 1/N*trapz([R R(1)]); % mean density of the solution

if (tpe == 1)
    eq4 = R_int - param; % find solution of a particular mean density
else
    eq4 = sqrt(2/N*trapz(([R R(1)] - R_int).^2)) - param; % find solution of a particular wave amplitude
end;

out = [eq1 eq2 eq3 eq4];

function out = Matrix(d,c,alpha,nu)
% matrix M-tilde defined in Eq. 57 of SI Appendix

rts = roots([1 -(d + (1+alpha^2)*(1-c))/(c) -d*(1+alpha^2)/(c)]);
R0 = min(rts);
C0 = (d - c*R0)/(1 + alpha^2);
out = [1/nu*(c - C0 - 1) -R0/nu 0; -1 -1/c -alpha/c; 0 alpha/c -1/c];

function out = func_dCalc(d,c,alpha,nu)
% function that finds the eigenvalue with the largest real part of the
% ... matrix M-tilde defined in Eq. 57 of SI Appendix, for given c and d.

M = Matrix(d,c,alpha,nu);
evals = eig(M);
out = max(real(evals));

function out = func(x,k,alpha,nu)
% function that finds the eigenvalue of the matrix M-tilde (defined in 
% ... Eq. 57 of SI Appendix) closest to 1i*k, for a given x = [c d]

c = x(1);
d = x(2);
M = Matrix(d,c,alpha,nu);
evals = eig(M);
out = [real(evals(2)); imag(evals(2))-k];


