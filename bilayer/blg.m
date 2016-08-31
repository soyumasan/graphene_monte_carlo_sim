% Monte Carlo simulation for spin current under a given Hamiltonian
% The boundary condition has been modified from what Privman uses
% It incorporates the magnetic field produced by electrons

clear
tstart = tic;

% Defining parameters
%--------------------------------------------------------------------------

dt = 0.005;                    % stepsize in time (fs)
Nsteps = 300000;              % Total no. of steps to be taken
effNsteps = 20000;            % No. of steps for which data is to be recorded
%Egrnd = 0.2*1.6e-19;         % first subband energy in J
%m = 0.050 * 9.1e-31;         % effective mass in kg
Nd = 1e18;                   % Donor ion conecntration in /m2 (assuming complete ionization)

%  Nd = 1e16;                 % for graphene


Lx = 5000;                    % channel length in nm
Ly = 10;                    % channel width in nm
Lz = 10;                     % Length along the confined direction in nm
N = round(Nd*Lx*Ly*1e-18);   % No. of electrons
kb = 1.38e-23;               % Boltzmann constant in SI units
T = 300;                     % Temperature in K


kT = kb*T;                   % Value of kT in SI units
phis = 0;                    % Source voltage in volts
%phid = 10.0;                 % Drain voltage in volts

phid = 4.0;

gridpoints = 500;            % No. of points in the grid for solving Poissons equation (Keep it even)
gridsize = Lx/gridpoints;    % grid spacing
I = 0;                       % Initialising the variable for current

dummy0 = 0;
dummy1 = 0;                  % dummy variables for random use
dummy2 = 0;
dummy3 = zeros(gridpoints,1);
gridE = zeros(1,gridpoints); % Allocating space for gridE vector

epsilonr = 2.3;            % relative permittivity of the material
g = -35 ;     			     % Lande g factor for electrons in the material concerned
mur = 1;                     % relative permeability of the material
factor1 = (gridsize * 1.6*100)/(epsilonr*8.85*Ly*Lz);       % prefactor used in Poisson equation solution
factor2 = g*4.615e-24;       % Factor for converting sigma to magnetic moment (2*gyromagnetic/hbar)

def_fac = 8.0;

acdefpot = def_fac*1.6e-19;      % Deformation potential for acoustic phonon scattering in SI units
hbar = 1.05e-34;             % Dirac's constant in SI units
longvel = 20000;              % Longitudinal velocity of sound in InGaAs m/sec 
%density = (5360+5667)/2;     % Density in SI units

density = 7.6138e-7;

nI = 6e13;                  % 3d Impurity density /m2 (dominated by dopant ions)

nI0 = 1e9;




% nI = 1e18;

D = 3.35e-1;                         % distance between graphene planes in nm;

nI1 = nI/(1e-9*D);                    % 3d impurity density in /m3

Nd1 = Nd/(1e-9*D);                   % 3d electron conc in /m3 



Z = 1;                       % Impurity charge

hbarwop = 0.196*1.6e-19;   % Optical phonon energy in SI units
DtK = 1.6e-19*def_fac*1e10;          % Optical coupling constant in SI units

eeta = 0.074*1.6e-29;        % eeta for rashba interaction in SI units
betadeff = eeta/5.3;         % Beta effective for dresselhaus interaction

Pscatt = 0;


Con_sp = 1;

%Con_sp = 0;

Con = 1;

%Con_sp = 0;

%Con_sp = 0;

%Con = 0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hbar = 1.05e-34;             % Dirac's constant in SI units

U1 = 0.5;                    % potential energy of the first layer in eV
U2 = 0.0;                    % potential energy of the second layer in eV
U = U1 - U2;


t = 2.7;                     % in eV
tperp = 0.35;                 % in eV 

%tperp = 3.5

vf = 1e6;                    % in m/sec                  

%U = 0.5;


D = 3.35e-1;    % distance between graphene planes in nano-metres


E1 = U/D;  % electric field in V/nanometers

U;

10*E1;    % in MV/cm

E6 = 10*E1;


zeta = 0.005; % (Rashba term from Physical Review B 80, 041405 (2009), by Ertler, Fabian et. al.) in meV/(V/nm)



del_1 = (1.6e-22)*2*zeta*E1/(hbar);




%A = (1.6e-22)*2*zeta*E/(hbar)

%B = 0;

%C = A;

Pscatt = 0;

m0 = 9.1e-31;

if (U < 0.14) 

  meff = (0.09*U + 0.043)*m0;

else

  meff = (1.6e-19)*tperp*(power((U*U + tperp*tperp), 1.5))/(2*U*(U*U + 2*tperp*tperp)*vf*vf);

end


%meff = meff * 1.6e-19;  % effective mass in Kg

%m

m = meff;

m11 = m/m0;

 
%kmin =  sqrt((U*U + 2*tperp*tperp)/(U*U + tperp*tperp))*(U/2*vf*hbar); % check the units
 
Egap = U*tperp/(sqrt(U*U + tperp*tperp)); % check the units,  I think it is in eV

%facEnergy = Egap/2 + (U1+U2)/2;

%facEnergy = 0;

%Energy = Egap/2 + (hbar*hbar/(2*m))*power( (k - kmin), 2) + (U1+U2)/2; % units should be in J


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


betakT = 1/(kT);

elcharge= 1.6e-19;

epsilon0 = 8.854e-12;    % in F/m

qs = 1e9;         % in 1/m

%qs = sqrt(4*pi*elcharge*elcharge*nI*betakT/epsilon0)

s = (betakT)*(hbar)*(hbar)*(qs*qs)/(2*m);

Ei1 = 0;

Ei2 = 0;

for i=1:1000

%Ei1 = 0.577 + (log(-s)) + (-s) + s*s/4

  Ei1 = Ei1 + exp(-s - 0.1*i)/(-s - 0.1*i);

end

J1 = -1 -(1+s)*exp(s)*Ei1; 


s = s/4;

for i=1:1000

%Ei1 = 0.577 + (log(-s)) + (-s) + s*s/4

  Ei2 = Ei2 + exp(-s - 0.1*i)/(-s - 0.1*i);

end



%Ei2 = 0.577 + real(log(-s)) + (-s) + s*s/4

  J2 = -1 -(1+s)*exp(s)*Ei2; 

delta = 0.00057e-3;     % in eV, check spin orbit coupling of graphene

gamma_el = 2*delta*(delta+2*Egap)/( (delta + Egap)*(2*delta + 3*Egap ) );

Prefac1 =  (elcharge^2*gamma_el/(epsilon0*Egap*elcharge))^2;

Prefac2 = sqrt(pi*kT/m);

%Prefac_e_e = 1.25*Nd1*Prefac1*Prefac2*J1; % both N and nI should be in 1/m^3;

Prefac_e_e = 1.25*Nd1*Prefac1*Prefac2*J1; % both N and nI should be in 1/m^3;

Prefac_e_imp = nI1*Prefac1*1.414*Prefac2*J2;


Prefac_e_e; 

Prefac_e_imp;


%spin_flip = Con_sp*(Prefac_e_e + Prefac_e_imp);

spin_flip = 0;

spin_flip;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Initialization
%--------------------------------------------------------------------------

x = Lx*(rand(1,N) - .5);                                     % x coordinate in nm
y = Ly*(rand(1,N) - .5);                                     % y coordinate in nm
dummy4 = realsqrt(-2*kT*reallog(rand(1,N))./m);              % velocity magnitude in m/sec (not carried forward)
dummy5 = 2*pi*rand(1,N);                                     % direction of velocity in rad (not carried forward)
Vx = dummy4.*cos(dummy5);                                    % x component of velocity in m/sec
Vy = dummy4.*sin(dummy5);                                    % y component of velocity in m/sec

dummy4 = acos(2*rand(1,N) - 1);                              % Choose theta
dummy5 = 2*pi*rand(1,N);                                     % Choose phi
Sz = cos(dummy4);                                            % Assign values to Sz, Sy and Sx
Sy = sin(dummy4).*sin(dummy5);                                  
Sx = sin(dummy4).*cos(dummy5);

gridx = -Lx/2 + gridsize/2 + (0:gridpoints-1)*gridsize;      % x coordinates of grid centres -- array

A = zeros(gridpoints,gridpoints);                            % the matrix for poisson's equation solution (later used for something else)
for dummy1 = 2:gridpoints-1
    A(dummy1,dummy1) = 2;
    A(dummy1,dummy1-1) = -1;
    A(dummy1,dummy1+1) = -1;
end
A(1,1) = 3;
A(1,2) = -1;
A(gridpoints,gridpoints) = 3;
A(gridpoints,gridpoints-1) = -1;
Ainv = inv(A);
                               
dummy3(1) = 2*phis;                                           % Poisson equation solving for the boundary part
dummy3(gridpoints) = 2*phid;
phibound = Ainv*dummy3;


% Acoustic phonon scattering rate in per sec

%acphonscat = (1.414*(m^(3/2))*kT*(acdefpot^2)*sqrt(Egrnd))/(pi*(hbar^4)*(longvel^2)*density);

acphonscat = (m*kT*(acdefpot^2))/((hbar^3)*(longvel^2)*density);

acphonscat


%% Optical phonon scattering rate in per sec

%Nop = 1/(exp(hbarwop/kT) - 1);
%optphonscatabs = ((DtK^2)*(m^1.5)*Nop*sqrt(Egrnd+hbarwop))/(1.414*pi*hbarwop*(hbar^2)*density);
%optphonscatem = ((DtK^2)*(m^1.5)*(Nop+1)*sqrt(Egrnd-hbarwop))/(1.414*pi*hbarwop*(hbar^2)*density);

%optphonscatabs = ((DtK^2)*(m)*Nop)/(2*hbarwop*(hbar)*density);


%%optphonscatem = ((DtK^2)*(m)*(Nop + 1))/(2*hbarwop*(hbar)*density);



%%optphonscatabs

%%optphonscatem




% Impurity scattering rate in per sec (CW formulation)

%b = (3/(4*pi*nI))^(1/3);

%b= (1/(pi*nI0))^(1/2);                  % in 2d

%nI


%epsb = ((1.6e-19)^2)/(2*epsilonr*8.85e-12*b);
%impscat = (Z^2)*pi*nI*(b^2)*sqrt(2*Egrnd/m);

%impscat = hbar*(Z^2)*pi*nI*(b^2)*(1/m);

%impscat 

%spin_flip



%%totalscat10 = acphonscat + optphonscatabs + optphonscatem;% + impscat
%%totalscat10

totalscat = acphonscat;%% + optphonscatabs + optphonscatem; % + impscat;% + spin_flip + 1e1;         % Total scattering rate
%%totalscat

%totalscat = 1e10;

%scatsteps = ceil((-1e15 * (log(rand(1,N))/totalscat))/dt);                 % Steps remaining till next scattering event


%scatratios = [acphonscat optphonscatabs optphonscatem impscat spin_flip]/totalscat;
%critnos = [scatratios(4) scatratios(3)+scatratios(4) scatratios(3)+scatratios(4)+scatratios(2) scatratios(3)+scatratios(4)+scatratios(2)+scatratios(1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


scatsteps = ceil((-1e15 * (log(rand(1,N))/totalscat))/dt);                 % Steps remaining till next scattering event


%scatratios[i] = [acphonscat[i] optphonscatabs[i] optphonscatem[i] impscat[i]]/totalscat[i];
%critnos_imp[i] = [scatratios(4) scatratios(3)+scatratios(4) scatratios(3)+scatratios(4)+scatratios(2) scatratios(3)+scatratios(4)+scatratios(2)+scatratios(1)];

%critnos_imp =  impscat/totalscat;

%%critnos_opt_em =  (optphonscatem)/totalscat;
%%critnos_opt_em

%%critnos_opt_abs =  (optphonscatabs + optphonscatem)/totalscat;
%%critnos_opt_abs

%%critnos_ac =  ( acphonscat + optphonscatabs + optphonscatem)/totalscat;
%%critnos_ac

%critnos_sp_fl =  ( acphonscat + optphonscatabs + optphonscatem + spin_flip )/totalscat;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%critnos

% Updatig grid and solving Poisson's equation
%--------------------------------------------------------------------------

gridindex = floor(x/gridsize + gridpoints/2 + 1);               % which grid the electron lies in

gridrho = (N/gridpoints)*ones(gridpoints,1) + calcrho(gridpoints,gridindex);



%for dummy2 = gridindex                                          % Calculation of charge density
%    gridrho(dummy2) = gridrho(dummy2) - 1;
%end

phicharge = factor1*Ainv*gridrho;                               % Poisson equation solving for charge distribution
phi = phicharge + phibound;                                     % Net voltage at the grid points

for dummy2 = 1:gridpoints-1                                     % Calculation of E field in kV/cm
    dummy3(dummy2) = (phi(dummy2+1) + phi(dummy2))/2;
end
dummy3(gridpoints) = phid;
gridE(2:gridpoints) = -1e4 * diff(dummy3)/gridsize;
gridE(1) = 1e4 * (phis - dummy3(1))/gridsize;

% Calculating spin distribution and magnetic field
%--------------------------------------------------------------------------

%gridSx = zeros(1,gridpoints);                                   % Allocating memory to variables
%gridSy = zeros(1,gridpoints);
%gridSz = zeros(1,gridpoints);

%for dummy2 = 1:N                                                % Calculation of grid Sx,Sy and Sz
%    dummy1 = gridindex(dummy2);                                 % Get the gridindex of the electron
%    gridSx(dummy1) = gridSx(dummy1) + Sx(dummy2);
 %   gridSy(dummy1) = gridSy(dummy1) + Sy(dummy2);
  %  gridSz(dummy1) = gridSz(dummy1) + Sz(dummy2);
%end

[gridSx gridSy gridSz] = calcgridS(gridpoints,gridindex,Sx,Sy,Sz);

gridN = N/gridpoints - gridrho';                                % Calculation of no. of electrons in each grid

% Magnetic field due to self interaction of the grids in SI units

gridBz =  -2.4*mur*1e-7*factor2.*gridSz.*(gridN.^2)./(Ly^3 * 1e-27);    
gridBy =  3*2.4*mur*1e-7*factor2.*gridSy.*(gridN.^2)./(Ly^3 * 1e-27);
gridBx = -2.4*mur*1e-7*factor2.*gridSx.*(gridN.^2)./(Ly^3 * 1e-27);

% Addition of magnetic field due to interaction between the grids in SI units

for dummy1 = 1:gridpoints
    dummy3 = -mur*2e-7*factor2*gridSz(dummy1)./(Ly.*(gridx(dummy1) - gridx + eps).^2 .* 1e-27);
    dummy3(dummy1) = 0;
    gridBz = gridBz + dummy3;
    
    dummy3 = -mur*2e-7*factor2*gridSx(dummy1)./(Ly.*(gridx(dummy1) - gridx + eps).^2 .* 1e-27);
    dummy3(dummy1) = 0;
    gridBx = gridBx - dummy3;
end
    
% Allocating memory to variables
%--------------------------------------------------------------------------

meanrho = zeros(gridpoints,1);
meanSx = zeros(1,gridpoints);
meanSy = zeros(1,gridpoints);
meanSz = zeros(1,gridpoints);
meanE = zeros(1,gridpoints);
meanBx = zeros(1,gridpoints);
meanBy = zeros(1,gridpoints);
meanBz = zeros(1,gridpoints);
meanPx = zeros(1,gridpoints);
meanPy = zeros(1,gridpoints);
meanPz = zeros(1,gridpoints);


meanPx_1 = zeros(1,gridpoints);
meanPy_1 = zeros(1,gridpoints);
meanPz_1 = zeros(1,gridpoints);
meanP_1 = zeros(1,gridpoints);





% Main simulation begins
%--------------------------------------------------------------------------


for dummy = 1:Nsteps
   
   if (mod(dummy,1000 ) == 0 )

    dummy

%  file_P2z = fopen('output_P2z.txt', 'w');
%  for z=1:gridpoints
%   Q(z) = realsqrt(Qz(z)*Qz(z) + Qy(z)*Qy(z) + Qx(z)*Qx(z) );
%  fprintf(file_P2z, '%8.6E  %8.6E %8.6E %8.6E %8.6E \n', z, Q(z)/num(z), Qz(z)/num(z), Qy(z)/num(z), Qx(z)/num(z) );
%  end

   end



 
    %disp(dummy)

% Flight calculations    

%dummy4 = (-1.6e-19/m)*(1e5*gridE(gridindex) + Vy.*gridBz(gridindex));   % Acceleration along x in m/sec2
%dummy5 = (1.6e-19/m)*(Vx.*gridBz(gridindex));                           % Acceleration along y in m/sec2
%x = x + 1e-6*dt*Vx + 0.5*1e-21*(dt^2)*dummy4 ;                          % in nm
%Vx = Vx + 1e-15*dt*dummy4;                                              % in m/sec
%y = y + 1e-6*dt*Vy + 0.5*1e-21*(dt^2)*dummy5;                           % in nm
%Vy = Vy + 1e-15*dt*dummy5;                                              % in m/sec



x = x + 1e-6*dt*Vx;                          % in nm

y = y + 1e-6*dt*Vy;                           % in nm



%dummy4 = (-1.6e-19/m)*(1e5*gridE(gridindex) + Vy.*gridBz(gridindex));   % Acceleration along x in m/sec2
%dummy5 = (1.6e-19/m)*(Vx.*gridBz(gridindex));                           % Acceleration along y in m/sec2

dummy4 = (-1.6e-19/m)*(1e5*gridE(gridindex));   % Acceleration along x in m/sec2

%dummy5 = (1.6e-19/m)*(Vx.*gridBz(gridindex));                           % Acceleration along y in m/sec2


Vx = Vx + 1e-15*dt*dummy4;                                              % in m/sec

Vy = Vy;                                              % in m/sec

%x = x + 1e-6*dt*Vx + 0.5*1e-21*(dt^2)*dummy4 ;                          % in nm

%y = y + 1e-6*dt*Vy + 0.5*1e-21*(dt^2)*dummy5;                           % in nm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for i=1:N

% cos_angle(i) = Vx(i)/sqrt(Vx(i)*Vx(i) + Vy(i)*Vy(i));

% sin_angle(i) = Vy(i)/sqrt(Vx(i)*Vx(i) + Vy(i)*Vy(i));

%%%%%%%%%%%%%%%%%%%%%%%%%% Spin updation parameters

 
%   B(i) = del_1*cos_angle(i);
%   A(i) = -del_1*sin_angle(i);
%   C(i) = del_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    
% Spin updation

%A = (2*m/(hbar^2))*(eeta*Vy - betadeff*(Vx - 1e-15*dt*0.5*dummy4));      % Calculating contributions from SO interaction in SI units
%B = (2*m/(hbar^2))*(betadeff*Vy - eeta*(Vx - 1e-15*dt*0.5*dummy4));      

%A = A - g*8.832e10*gridBx(gridindex);                           % Defining parameters for updating spin
%B = B - g*8.832e10*gridBy(gridindex);
%C = -g*8.832e10*gridBz(gridindex);




dummy4 = Sz;                                                    % The differential way of spin updating
dummy5 = Sx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  C =  (1.6e-22)*2*zeta*E/(hbar);

%  A =  (1.6e-22)*2*zeta*E/(hbar);

%  C = 0;
 
%  A = 0; 

%  B = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Sz = Sz + 1e-15*dt*(A.*Sy - B.*Sx);
%Sx = Sx + 1e-15*dt*(B.*dummy4 - C.*Sy);
%Sy = Sy + 1e-15*dt*(C.*dummy5 - A.*dummy4);



for i = 1:N

%i
%A(i)
%B(i)
%C(i)

%Sz(i) = Sz(i) + 1e-15*dt*(A(i)*Sy(i) - B(i)*Sx(i));
%Sx(i) = Sx(i) + 1e-15*dt*(B(i)*dummy4(i) - C(i)*Sy(i));
%Sy(i) = Sy(i) + 1e-15*dt*(C(i)*dummy5(i) - A(i)*dummy4(i));



  cos_angle(i) = Vx(i)/sqrt(Vx(i)*Vx(i) + Vy(i)*Vy(i));

  sin_angle(i) = Vy(i)/sqrt(Vx(i)*Vx(i) + Vy(i)*Vy(i));

%%%%%%%%%%%%%%%%%%%%%%%%%% Spin updation parameters

 
   B(i) = del_1*cos_angle(i);
   A(i) = -del_1*sin_angle(i);
   C(i) = del_1;



Sz_ori(i) = Sz(i);

Sx_ori(i) = Sx(i);

Sy_ori(i) = Sy(i);


Sz(i) = Sz(i) + 1e-15*dt*(A(i)*Sy_ori(i) - B(i)*Sx_ori(i));

Sx(i) = Sx(i) + 1e-15*dt*(B(i)*Sz_ori(i));

Sy(i) = Sy(i) + 1e-15*dt*(- A(i)*Sz_ori(i));


end






dummy4 = sqrt(Sx.^2 + Sy.^2 + Sz.^2);
Sx = Sx./dummy4;
Sy = Sy./dummy4;
Sz = Sz./dummy4;

% Imposing boundary conditions (assuming specular walls along y)
%--------------------------------------------------------------------------

indices = find(abs(y)>Ly/2);                                    % indices of electrons with |y| > 500
y(indices) = (2*(y(indices) > 0) - 1).*Ly - y(indices);         % updating y and Vy acc. to specular reflection
Vy(indices) = (2*(y(indices) < 0) - 1).*abs(Vy(indices));    


indices = find(x<-Lx/2);                                        % indices of electrons with x < 500
x(indices) = x(indices) + Lx;                                   % Introduce new electrons in accordance with PBC
dummy6 = -2*kT*log(rand(1,length(indices)))./m;                 % magnitude of velocity square
Vy(indices) = sqrt(dummy6).*sin(pi*rand(1,length(indices)) - pi/2);
Vx(indices) = -sqrt(dummy6 - Vy(indices).^2);
dummy6 = acos(2*rand(1,length(indices)) - 1);                   % Choose theta
dummy7 = 2*pi*rand(1,length(indices));                          % Choose phi
Sz(indices) = cos(dummy6);                                      % Assign values to Sz, Sy and Sx
Sy(indices) = sin(dummy6).*sin(dummy7);                                  
Sx(indices) = sin(dummy6).*cos(dummy7);

if dummy>Nsteps - effNsteps                                     % Calculation of no. of electrons that are injected by the drain
    I = I - length(indices);
end

indices = find(x>Lx/2);                                         % indices of electrons with x > 500
x(indices) = x(indices) - Lx;                                   % Introduce new electrons in accordance with PBC


dummy6 = -2*kT*log(rand(1,length(indices)))./m;                 % magnitude of velocity square
Vy(indices) = sqrt(dummy6).*sin(pi*rand(1,length(indices)) - pi/2);
Vx(indices) = sqrt(dummy6 - Vy(indices).^2);
Sz(indices) = 1;                                                % updating the spin of the electrons introduced
Sx(indices) = 0;
Sy(indices) = 0;

if dummy>Nsteps - effNsteps                                     % Calulation of no. of electrons that enter the drain
    I = I + length(indices);
end

    
% Updatig grid and solving Poisson's equation
%--------------------------------------------------------------------------

gridindex = floor(x/gridsize + gridpoints/2 + 1);               % which grid the electron lies in

gridrho = (N/gridpoints)*ones(gridpoints,1) + calcrho(gridpoints,gridindex);


%for dummy2 = gridindex                                          % Calculation of charge density
%    gridrho(dummy2) = gridrho(dummy2) - 1;
%end

phicharge = factor1*Ainv*gridrho;                               % Poisson equation solving for charge distribution
phi = phicharge + phibound;                                     % Net voltage at the grid points

for dummy2 = 1:gridpoints-1                                     % Calculation of E field in kV/cm
    dummy3(dummy2) = (phi(dummy2+1) + phi(dummy2))/2;
end
dummy3(gridpoints) = phid;
gridE(2:gridpoints) = -1e4 * diff(dummy3)/gridsize;
gridE(1) = 1e4 * (phis - dummy3(1))/gridsize;


% Calculating spin distribution and magnetic field
%--------------------------------------------------------------------------

%gridSx = zeros(1,gridpoints);                                   % Allocating memory to variables
%gridSy = zeros(1,gridpoints);
%gridSz = zeros(1,gridpoints);

%for dummy2 = 1:N                                                % Calculation of grid Sx,Sy and Sz
%    dummy1 = gridindex(dummy2);                                 % Get the gridindex of the electron
%    gridSx(dummy1) = gridSx(dummy1) + Sx(dummy2);
 %   gridSy(dummy1) = gridSy(dummy1) + Sy(dummy2);
  %  gridSz(dummy1) = gridSz(dummy1) + Sz(dummy2);
%end

[gridSx gridSy gridSz] = calcgridS(gridpoints,gridindex,Sx,Sy,Sz);


gridN = N/gridpoints - gridrho';                                % Calculation of no. of electrons in each grid

% Magnetic field due to self interaction of the grids in SI units

gridBz =  -2.4*mur*1e-7*factor2.*gridSz.*(gridN.^2)./(Ly^3 * 1e-27);    
gridBy =  3*2.4*mur*1e-7*factor2.*gridSy.*(gridN.^2)./(Ly^3 * 1e-27);
gridBx = -2.4*mur*1e-7*factor2.*gridSx.*(gridN.^2)./(Ly^3 * 1e-27);

% Addition of magnetic field due to interaction between the grids in SI units

for dummy1 = 1:gridpoints
    dummy3 = -mur*2e-7*factor2*gridSz(dummy1)./(Ly.*(gridx(dummy1) - gridx + eps).^2 .* 1e-27);
    dummy3(dummy1) = 0;
    gridBz = gridBz + dummy3;
    
    dummy3 = -mur*2e-7*factor2*gridSx(dummy1)./(Ly.*(gridx(dummy1) - gridx + eps).^2 .* 1e-27);
    dummy3(dummy1) = 0;
    gridBx = gridBx - dummy3;
end


% Scattering calculations


indices = find(scatsteps==1);                                   % Indices of the electrons undergoing scattering

scatsteps = scatsteps - 1;                                      % Update scatsteps
scatsteps(indices) = ceil((-1e15 * (reallog(rand(1,length(indices)))/totalscat))/dt); 

for dummy2 = indices
    dummy0 = rand;
%    if dummy0 < critnos(1)                                      % Impurity scattering
 
      
     %if dummy0 < critnos_imp
       
      %    Pscatt = 1;


      %  dummy1 = atan2(Vy(dummy2),Vx(dummy2));                  % Calculating the current angle of velocity
      %  Vx(dummy2) = sqrt(Vx(dummy2)^2 + Vy(dummy2)^2);         % Magnitude of velocity
      %  dummy0 = rand;
      %  Vy(dummy2) = acos((((Egrnd + 0.5*m*Vx(dummy2)^2/epsb)^2)*dummy0 - 1)/(((Egrnd + 0.5*m*Vx(dummy2)^2/epsb)^2)*dummy0 + 1)) * (2*(rand>0.5) - 1);       % thetar
      %  dummy1 = dummy1 + Vy(dummy2);                           % Final angle of velocity
      %  Vy(dummy2) = Vx(dummy2) * sin(dummy1);                  % updating Vy
      %  Vx(dummy2) = Vx(dummy2) * cos(dummy1);                  % updating Vx
        

       % Vx(dummy2) = sqrt(Vx(dummy2)^2 + Vy(dummy2)^2);                     % store the magnitude of V in Vx
     %   Vy(dummy2) = Vx(dummy2)*sin(2*pi*rand);                             % update Vy
      %  Vx(dummy2) = sqrt(Vx(dummy2)^2 - Vy(dummy2)^2)*(2*(rand>0.5) - 1);  % update Vx
        



%    elseif ( (dummy0 < critnos_opt_em ) && ( (Vx(dummy2)^2 + Vy(dummy2)^2 - 2*hbarwop/m) > 0) )    % Optical phonon emmission
        
     %%if (dummy0 < critnos_opt_em ) && (Vx(dummy2)^2 + Vy(dummy2)^2 - 2*hbarwop/m > 0)     % Optical phonon emmission
 

       %% Pscatt = 2;


        %%Vx(dummy2) = sqrt(Vx(dummy2)^2 + Vy(dummy2)^2 - 2*hbarwop/m);       % store the updated magnitude of v in Vx
        %%Vy(dummy2) = Vx(dummy2)*sin(2*pi*rand);                             % update Vy
        %%Vx(dummy2) = sqrt(Vx(dummy2)^2 - Vy(dummy2)^2)*(2*(rand>0.5) - 1);  % update Vx
        
%    elseif dummy0 < critnos(3)                                  % Optical phonon absorption
 
       %%elseif (dummy0 < critnos_opt_abs)                                  % Optical phonon absorption
      
        %% Pscatt = 3;

 
       %% Vx(dummy2) = sqrt(Vx(dummy2)^2 + Vy(dummy2)^2 + 2*hbarwop/m);       % store the updated magnitude of v in Vx
       %% Vy(dummy2) = Vx(dummy2)*sin(2*pi*rand);                             % update Vy
       %% Vx(dummy2) = sqrt(Vx(dummy2)^2 - Vy(dummy2)^2)*(2*(rand>0.5) - 1);  % update Vx
        
 %   elseif dummy0 < critnos(4)                                                        % Acoustic phonon scattering
  

      if (dummy0 < 1 )                   %%critnos_ac for this case is 1                                     % Acoustic phonon scattering    
      
        Pscatt = 4;


        Vx(dummy2) = sqrt(Vx(dummy2)^2 + Vy(dummy2)^2);                     % store the magnitude of V in Vx
        Vy(dummy2) = Vx(dummy2)*sin(2*pi*rand);                             % update Vy
        Vx(dummy2) = sqrt(Vx(dummy2)^2 - Vy(dummy2)^2)*(2*(rand>0.5) - 1);  % update Vx
        

      %elseif (dummy0 < critnos_sp_fl )    % spin flip

       %Pscatt = 5;

       %Sz(dummy2) = -Sz(dummy2);
       %Sy(dummy2) = -Sy(dummy2);
       %Sx(dummy2) = -Sx(dummy2);



     else

        Pscatt = 6;

       Vx(dummy2) = Vx(dummy2);
       Vy(dummy2) = Vy(dummy2);

  


    end
    
    
end

if dummy > (Nsteps-effNsteps)
    
    meanrho = meanrho + gridrho;
    meanSx = meanSx + gridSx;
    meanSy = meanSy + gridSy;
    meanSz = meanSz + gridSz;
    meanBx = meanBx + gridBx;
    meanBy = meanBy + gridBy;
    meanBz = meanBz + gridBz;
    meanPx = meanPx + gridSx./gridN;
    meanPy = meanPy + gridSy./gridN;
    meanPz = meanPz + gridSz./gridN;
    


meanrho_1 = meanrho/(dummy-Nsteps+effNsteps);
meanSx_1 = meanSx/(dummy-Nsteps+effNsteps);
meanSy_1 = meanSy/(dummy-Nsteps+effNsteps);
meanSz_1 = meanSz/(dummy-Nsteps+effNsteps);

meanN_1 = N/gridpoints - meanrho_1';
meanPx_1 = meanPx/(dummy-Nsteps+effNsteps);
meanPy_1 = meanPy/(dummy-Nsteps+effNsteps);
meanPz_1 = meanPz/(dummy-Nsteps+effNsteps);
meanBx_1 = meanBx/(dummy-Nsteps+effNsteps);
meanBy_1 = meanBy/(dummy-Nsteps+effNsteps);
meanBz_1 = meanBz/(dummy-Nsteps+effNsteps);


%end

meanP_1 = realsqrt(meanPx_1.^2 + meanPy_1.^2 + meanPz_1.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ((mod(dummy,1000 ) < 0.00001 ))


file_P3z = fopen('output_P2z_2.txt', 'w');

for z=1:gridpoints

%   Q(z) = realsqrt(Qz(z)*Qz(z) + Qy(z)*Qy(z) + Qx(z)*Qx(z) );

%  fprintf(file_P2z, '%8.6E  %8.6E %8.6E %8.6E %8.6E \n', z, Q(z)/num(z), Qz(z)/num(z), Qy(z)/num(z), Qx(z)/num(z) );


%  fprintf(file_P2z, '%8.6E  %8.6E \n', z, Q(z)/num(z) );


   fprintf(file_P3z, '%8.6E  %8.6E %8.6E %8.6E %8.6E  \n', z, meanP_1(z), meanPx_1(z), meanPy_1(z), meanPz_1(z));


end

fclose(file_P3z);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end



end




end

meanrho = meanrho/effNsteps;
meanSx = meanSx/effNsteps;
meanSy = meanSy/effNsteps;
meanSz = meanSz/effNsteps;

meanN = N/gridpoints - meanrho';
meanPx = meanPx/effNsteps;
meanPy = meanPy/effNsteps;
meanPz = meanPz/effNsteps;
meanBx = meanBx/effNsteps;
meanBy = meanBy/effNsteps;
meanBz = meanBz/effNsteps;

meanP = realsqrt(meanPx.^2 + meanPy.^2 + meanPz.^2);

meanphi = factor1*Ainv*meanrho + phibound;                      % Net voltage at the grid points

for dummy2 = 1:gridpoints-1                                     % Calculation of E field in kV/cm
    dummy3(dummy2) = (meanphi(dummy2+1) + meanphi(dummy2))/2;
end
dummy3(gridpoints) = phid;
meanE(2:gridpoints) = -1e4 * diff(dummy3)/gridsize;
meanE(1) = 1e4 * (phis - dummy3(1))/gridsize;

I = (0.5*1.6e-19*I)/(effNsteps*dt*1e-15);                       % Calculation of current

% Calculation of spin relxation length

indices = find(abs(meanP-exp(-1))<0.02);
Lsabs = mean(gridx(indices)) + Lx/2;

xfit = (gridx + Lx/2);
coefficients = polyfit(xfit,log(meanP),1);
Lsfit = -1/coefficients(1);
S0 = exp(coefficients(2));

% Fitting enforcing the condition of S0 = 1

Lsdiff = -1/mean(diff(log(meanP))./diff(gridx));

file_P2z = fopen('output_P2z_1.txt', 'w');


%num

for z=1:gridpoints

%   Q(z) = realsqrt(Qz(z)*Qz(z) + Qy(z)*Qy(z) + Qx(z)*Qx(z) );

%  fprintf(file_P2z, '%8.6E  %8.6E %8.6E %8.6E %8.6E \n', z, Q(z)/num(z), Qz(z)/num(z), Qy(z)/num(z), Qx(z)/num(z) );


%  fprintf(file_P2z, '%8.6E  %8.6E \n', z, Q(z)/num(z) );


   fprintf(file_P2z, '%8.6E  %8.6E %8.6E %8.6E %8.6E  \n', z, meanP(z), meanPx(z), meanPy(z), meanPz(z));


end

fclose(file_P2z);





timetaken = toc(tstart);
timetaken = timetaken/60;
disp(timetaken)

%cd data
save 1e18g-35
%cd ..