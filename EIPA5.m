%PA5 - winter 2022 ELEC 4700 Harmonic Wave Equation in 2D FD and Modes
%11 feb 2022
%Alina Jacobson

clear all;
close all;

set(0,'DefaultFigureWindowStyle','docked')

global C

%initialize global constants
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                 % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      %metres (32.1740 ft) per s²

% i - Initially set nx = ny = 50;
nx=50;
ny=50;

% nx=65;          %with different values - larger enteries 3000
% ny=50;          %and non-zeros increases 15120

% i -Create a sparse G matrix (nx*ny,nx*ny) in size.
G=sparse(nx*ny,nx*ny);
newG = zeros(nx,ny);       %for remapping - initial matrix all zeros

%loop through to set the BC’s nodes diagonal value G(m,m) to 1 
%and all other entries to 0 (G(m,:))
for i=1:nx      %nx

    for j=1:ny  %ny
        
%         (i,j) -> n
%         (i-1,j) -> nxm
%         (i+1,j) -> nxp
%         (i,j-1) -> nym
%         (i,j+1) -> nyp

        %mapping equation n=j+(i-1)*ny
        n = j + (i-1) * ny;     
        nxm = j + (i-2) * ny;	
        nxp = j + i * ny;        
        nym = j-1 + (i-1) * ny; 
        nyp = j+1 + (i-1) * ny; 

        % ii -set BC’s nodes diagonal value G(m,m) to 1 else other entries to 0 (G(m,:))
        if (i==1 || i==nx )
            G(n,n)=1; 
        elseif(j==1||j==ny)
            G(n,n)=1;
        else
            
        % iii - For the bulk nodes the G entries are set by the FD equation
        %previous entry for nx and ny
        
        %using -4 in the G matrix diagonal - inner region equation
        G(n,n) = -4;
     
        G(n,nxm) = 1;
        G(n,nxp) = 1; 
        G(n,nym) = 1;
        G(n,nyp) = 1;
        
        
        % ix change to -2
        if i > 10 && i < 20 && j > 10 && j < 20
            G(n,n) = -2;
        end
      
        end
    end
end

% iV - plot G matrix - square matrix - number of unknowns by number of unknows
% spy() plots all the non-zeros ( zoomed in shows 5 enteries)
figure(1)
spy(G)                          % non-zeros=11520   number of enteries = 2500
title('Matrix G - using spy(G)')

% v - Use [E,D] = eigs(g,9,’SM’) to get 9 eigenvectors and values. 
% E is a matrix of vectors and D’s diagonal is the eigenvalues.
%plot the  Eigen Values
[E,D] = eigs(G, 9, 'SM');

figure(2)
ev=diag(D);          %diagonal is the eigenvalues
plot(ev, '*');
title('Matrix G - Eigen values using eigs')


%vii. Plot the eigenvectors using surf(). 
%will need to remap them onto a matrix that is (nx,ny) in size. 
for z = 1:9                 %outputs 3x3
    matrixE=E(:,z);         %define matrix E with value range
    
    %loop again for nx
    for i=1:nx
        
        %loop again for ny
        for j=1:ny
            
            %remap using same map equation as above
            n = j + (i-1) * ny;  
            newG(i,j)= matrixE(n);     %newG matrix is given matrixE defined above
        end
        
        
    end
    
    %output plot on same page 3x3 with 9 surface plots
    plotALL = z;
    figure(3)
    subplot(3,3,plotALL);
    surf(newG);
    title('Eigenvectors using surf()')    
    
end

        
