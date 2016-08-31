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



Lx = 5000;                    % channel length in nm
Ly = 10;                    % channel width in nm
Lz = 10;                     % Length along the confined direction in nm
N = round(Nd*Lx*Ly*1e-18);   % No. of electrons
kb = 1.38e-23;               % Boltzmann constant in SI units
T = 300;                     % Temperature in K
kT = kb*T;                   % Value of kT in SI units
phis = 0;                    % Source voltage in volts
phid = 4.0;                 % Drain voltage in volts
gridpoints = 500;            % No. of points in the grid for solving Poissons equation (Keep it even)
gridsize = Lx/gridpoints;    % grid spacing
I = 0;                       % Initialising the variable for current

%N = 20;

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

acdefpot = 8.0*1.6e-19;      % Deformation potential for acoustic phonon scattering in SI units
hbar = 1.05e-34;             % Dirac's constant in SI units
longvel = 20000;              % Longitudinal velocity of sound in InGaAs m/sec 
%density = (5360+5667)/2;     % Density in SI units

density = 7.6138e-7;

nI = 6e13;                  % 3d Impurity density /m3 (dominated by dopant ions)
Z = 1;                       % Impurity charge

hbarwop = 0.196*1.6e-19;   % Optical phonon energy in SI units
DtK = 1.6e-19*8e10;          % Optical coupling constant in SI units

eeta = 0.074*1.6e-29;        % eeta for rashba interaction in SI units
betadeff = eeta/5.3;         % Beta effective for dresselhaus interaction

Pscatt = 0;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hbar = 1.05e-34;             % Dirac's constant in SI units

U1 = 0.0;
U2 = 0.0;                    % potential energy of the second layer in eV
U = U1 - U2;


t = 2.7;                     % in eV
tperp = 3.5;                 % in eV 

vf = 1e6;                    % in m/sec                  

D = 3.35e-1;    % distance between graphene planes in nano-metres


E = U/D;  % electric field in V/nanometers

zeta = 0.005; % (Rashba term from Physical Review B 80, 041405 (2009), by Ertler, Fabian et. al.) in meV/(V/nm)


Pscatt = 0;

m0 = 9.1e-31;

if (U < 0.14) 

  meff = (0.09*U + 0.043)*m0;

else

  meff = (1.6e-19)*tperp*(power((U*U + tperp*tperp), 1.5))/(2*U*(U*U + 2*tperp*tperp)*vf*vf);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%Initialization
%--------------------------------------------------------------------------

x = Lx*(rand(1,N) - .5);                                     % x coordinate in nm
y = Ly*(rand(1,N) - .5);                                     % y coordinate in nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:N

x(i) = Lx*(rand - 0.5);

y(i) = Lx*(rand - 0.5);


mag(i) = -kT*reallog(rand)/vf;
angle(i) = 2*pi*rand;

en(i) = vf*mag(i);                                           % energy

kx(i) = mag(i)*cos(angle(i))/hbar;                                    % x component of wavevector in m/sec
ky(i) = mag(i)*sin(angle(i))/hbar;                                    % y component of wavevector in m/sec

modk(i) = sqrt(kx(i)*kx(i) + ky(i)*ky(i));


Vx(i) = vf*cos(angle(i));
Vy(i) = vf*sin(angle(i));


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dummy4 = acos(2*rand(1,N) - 1);                              % Choose theta
dummy5 = 2*pi*rand(1,N);                                     % Choose phi
Sz = cos(dummy4);                                            % Assign values to Sz, Sy and Sx
Sy = sin(dummy4).*sin(dummy5);                                  
Sx = sin(dummy4).*cos(dummy5);

gridx = -Lx/2 + gridsize/2 + (0:gridpoints-1)*gridsize;      % x coordinates of grid centres

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



%%%%%%%%%%%%%
acphonscatinitop = fopen('acphonscatinitop1.txt', 'a');
%%optphonscateminitop =fopen('optphonscateminitop1.txt', 'a');
%%optphonscatabsinitop =fopen('optphonscateminitop1.txt', 'a');


%%%%%%%%%%%%%





for i=1:N

en(i) = hbar*vf*sqrt(kx(i)*kx(i) + ky(i)*ky(i));

% Acoustic phonon scattering rate in per sec

acphonscat(i) = (en(i)/(vf*vf))*(kT*(acdefpot^2))/((hbar^3)*(longvel^2)*density);
fprintf(acphonscatinitop, '%8.6E  \n', acphonscat(i));


% Optical phonon scattering rate in per sec

%%Nop = 1/(exp(hbarwop/kT) - 1);
%%optphonscatabs(i) = ((DtK^2)*Nop*(en(i)+hbarwop))/(2*hbarwop*(hbar)*density*vf*vf);
%%fprintf(optphonscatabsinitop, '%8.6E  \n', optphonscatabs(i));

%%if ( (en(i)-hbarwop) > 0.00000001 )

%%optphonscatem(i) = ((DtK^2)*(Nop+1)*(en(i)-hbarwop))/(2*hbarwop*(hbar)*density*vf*vf);
%%fprintf(optphonscateminitop, '%8.6E  \n', optphonscatem(i));

  
%%else  

%%optphonscatem(i) = 0;

end


% Impurity scattering rate in per sec (CW formulation)

b = (3/(4*pi*nI))^(1/3);
epsb = ((1.6e-19)^2)/(2*epsilonr*8.85e-12*b);


totalscat(i) = acphonscat(i);%% + optphonscatabs(i) + optphonscatem(i); %%+ impscat(i) + 1e1;         % Total scattering rate

scatsteps(i) = ceil((-1e15 * (log(rand)/totalscat(i)))/dt);                 % Steps remaining till next scattering event


%%critnos_opt_em(i) =  (optphonscatem(i))/totalscat(i);

%%critnos_opt_abs(i) =  (optphonscatabs(i) + optphonscatem(i))/totalscat(i);

critnos_ac(i) = 1;%% ( acphonscat(i) + optphonscatabs(i) + optphonscatem(i))/totalscat(i);





end


% Updatig grid and solving Poisson's equation
%--------------------------------------------------------------------------

gridindex = floor(x/gridsize + gridpoints/2 + 1);               % which grid the electron lies in

gridrho = (N/gridpoints)*ones(gridpoints,1) + calcrho(gridpoints,gridindex);


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




%%%%%%%%%%%%%
acphonscatop = fopen('acphonscatop1.txt', 'a');
%%optphonscatemop =fopen('optphonscatemop1.txt', 'a');
%%optphonscatabsop =fopen('optphonscatemop1.txt', 'a');


acphonscatoplast = fopen('acphonscatop2.txt', 'w');
%%optphonscatemoplast =fopen('optphonscatemop2.txt', 'w');
%%optphonscatabsoplast =fopen('optphonscatemop2.txt', 'w');
%%%%%%%%%%%%%


% Main simulation begins
%--------------------------------------------------------------------------


for dummy = 1:Nsteps
   
      if (mod(dummy,1000 ) == 0 )

    dummy  %%%%%%%%%%%%status tracker


   end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  for i=1:N

kx(i) = kx(i) - 1.6e-19*(1e-15*dt)*(1e5*gridE(real(gridindex(i))) )/(hbar);

ky(i) = ky(i);


Vx(i) = vf*(kx(i))/(sqrt(kx(i)*kx(i) + ky(i)*ky(i)));  

qw =  (Vx(i))^2;

Vy(i) = sqrt(vf^2 - qw);


x(i) = x(i) + 1e-6*dt*Vx(i) ;

y(i) = y(i) + 1e-6*dt*Vy(i);                                                      % in nm



  end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





    
% Spin updation

dummy4 = Sz;                                                    % The differential way of spin updating
dummy5 = Sx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  C =  (1.6e-22)*2*zeta*E/(hbar);

  A =  (1.6e-22)*2*zeta*E/(hbar);

  B = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Sz = Sz + 1e-15*dt*(A.*Sy - B.*Sx);
Sx = Sx + 1e-15*dt*(B.*dummy4 - C.*Sy);
Sy = Sy + 1e-15*dt*(C.*dummy5 - A.*dummy4);

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


for i=indices

x(i) = x(i) + Lx;                                   % Introduce new electrons in accordance with PBC

mag(i) = -kT*reallog(rand)/vf;
angle(i) = 2*pi*rand;


kx(i) = mag(i)*cos(angle(i))/hbar;                                    % x component of wavevector in m/sec
ky(i) = mag(i)*sin(angle(i))/hbar;                                    % y component of wavevector in m/sec

Vx(i) = -vf*abs(cos(angle(i)));
Vy(i) = vf*sin(angle(i));


end


dummy6 = acos(2*rand(1,length(indices)) - 1);                   % Choose theta
dummy7 = 2*pi*rand(1,length(indices));                          % Choose phi
Sz(indices) = cos(dummy6);                                      % Assign values to Sz, Sy and Sx
Sy(indices) = sin(dummy6).*sin(dummy7);                                  
Sx(indices) = sin(dummy6).*cos(dummy7);

if dummy>Nsteps - effNsteps                                     % Calculation of no. of electrons that are injected by the drain
    I = I - length(indices);
end

indices = find(x>Lx/2);                                         % indices of electrons with x > 500


for i=indices

x(i) = x(i) - Lx;                                   % Introduce new electrons in accordance with PBC

mag(i) = -kT*reallog(rand)/vf;
angle(i) = 2*pi*rand;

kx(i) = mag(i)*cos(angle(i))/hbar;                                    % x component of wavevector in m/sec
ky(i) = mag(i)*sin(angle(i))/hbar;                                    % y component of wavevector in m/sec

Vx(i) = vf*abs(cos(angle(i)));
Vy(i) = vf*sin(angle(i));


end




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





for i=indices

en(i) = hbar*vf*sqrt(kx(i)*kx(i) + ky(i)*ky(i));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acphonscat(i) = (en(i)/(vf*vf))*(kT*(acdefpot^2))/((hbar^3)*(longvel^2)*density);


%%%%%
    if (dummy>(Nsteps-effNsteps))
        fprintf(acphonscatop, '%8.6E  %8.6E  \n',i,acphonscat(i));
        fprintf(acphonscatoplast, '%8.6E  %8.6E  \n',i,acphonscat(i));
    end

%%%%%


% Optical phonon scattering rate in per sec

%%Nop = 1/(exp(hbarwop/kT) - 1);
%%optphonscatabs(i) = ((DtK^2)*Nop*(en(i)+hbarwop))/(2*hbarwop*(hbar)*density*vf*vf);


%%%%%
   %% if (dummy>(Nsteps-effNsteps))
       %% fprintf(optphonscatabsop, '%8.6E  %8.6E  \n',i,optphonscatabs(i));
       %% fprintf(optphonscatabsoplast, '%8.6E  %8.6E  \n',i,optphonscatabs(i));
    %%end

%%%%%


%%if ( (en(i)-hbarwop) > 0.00000001 ) 

%%optphonscatem(i) = ((DtK^2)*(Nop+1)*(en(i)-hbarwop))/(2*hbarwop*(hbar)*density*vf*vf);

 

%% else 

 %%optphonscatem(i) = 0;

%%end

%%%%%
   %% if (dummy>(Nsteps-effNsteps))
    %%    fprintf(optphonscatemop, '%8.6E  %8.6E  \n',i,optphonscatem(i));
    %%    fprintf(optphonscatemoplast, '%8.6E  %8.6E  \n',i,optphonscatem(i));
   %% end

%%%%%

% Impurity scattering rate in per sec (CW formulation)

b = (3/(4*pi*nI))^(1/3);
epsb = ((1.6e-19)^2)/(2*epsilonr*8.85e-12*b);


totalscat(i) = acphonscat(i);%% + optphonscatabs(i) + optphonscatem(i); %% + impscat(i) + 1e1;         % Total scattering rate




scatsteps(i) = ceil((-1e15 * (log(rand)/totalscat(i)))/dt);                 % Steps remaining till next scattering event

%%critnos_opt_em(i) =  (optphonscatem(i))/totalscat(i);

%%critnos_opt_abs(i) =  (optphonscatabs(i) + optphonscatem(i))/totalscat(i);

critnos_ac(i) =  1;%%( acphonscat(i) + optphonscatabs(i) + optphonscatem(i))/totalscat(i);





end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














%scatsteps(i) = ceil((-1e15 * (reallog(rand(1,length(indices)))/totalscat(i)))/dt); 





for dummy2 = indices
    dummy0 = rand;
           
    %%if (   (dummy0 < critnos_opt_em(dummy2)) &&  (sqrt(kx(dummy2)^2 + ky(dummy2)^2) > hbarwop/vf ) )    % Optical phonon emmission


     %%    Pscatt = 2;



   %%    modkprime(dummy2) = sqrt(kx(dummy2)^2  + ky(dummy2)^2)  - hbarwop/vf;
    %%   ky(dummy2) = modkprime(dummy2)*sin(2*pi*rand);
     %%  kx(dummy2) =  sqrt(modkprime(dummy2)^2 - ky(dummy2)^2)*(2*(rand>0.5) - 1);




        
   %% elseif dummy0 < critnos_opt_abs(dummy2)                                  % Optical phonon absorption
       
       %%  Pscatt = 3;

 
      %% modkprime(dummy2) = sqrt(kx(dummy2)^2  + ky(dummy2)^2)  + hbarwop/vf;
      %% ky(dummy2) = modkprime(dummy2)*sin(2*pi*rand);
       %%kx(dummy2) =  sqrt(modkprime(dummy2)^2 - ky(dummy2)^2)*(2*(rand>0.5) - 1);





    if dummy0 < 1                                                        % Acoustic phonon scattering
        
        Pscatt = 4;




       modkprime(dummy2) = sqrt(kx(dummy2)^2  + ky(dummy2)^2);
       ky(dummy2) = modkprime(dummy2)*sin(2*pi*rand);
       kx(dummy2) =  sqrt(modkprime(dummy2)^2 - ky(dummy2)^2)*(2*(rand>0.5) - 1);


        
     else

        Pscatt = 5;

       kx(dummy2) = kx(dummy2);
       ky(dummy2) = ky(dummy2);

  


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




file_P3z = fopen('output_P2z_2.txt', 'w');

for z=1:gridpoints


   fprintf(file_P3z, '%8.6E  %8.6E %8.6E %8.6E %8.6E  \n', z, meanP_1(z), meanPx_1(z), meanPy_1(z), meanPz_1(z));


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

   fprintf(file_P2z, '%8.6E  %8.6E %8.6E %8.6E %8.6E  \n', z, meanP(z), meanPx(z), meanPy(z), meanPz(z));


end

fclose(file_P2z);





timetaken = toc(tstart);
timetaken = timetaken/60;
disp(timetaken)

%cd data
save 1e18g-35
%cd ..
