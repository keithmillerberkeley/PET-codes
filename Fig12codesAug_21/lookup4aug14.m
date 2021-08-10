
          %  file lookup4aug14.m

   %  called from randPET3D.m
      %  For a single beam, with angle phi we compute the kernel function      
      %  K(syy1,syy2,phi) on the screen coords syy1,syy2
      %  by table lookup from the precomputed TBLKPHI
      % I passed syy1 syy2 as globals

    function Hts = lookup4aug14(phi)  % I passed syy1 syy2 as globals


% phi is scalar;  syy1,syy2 are the screen coords of (recon pts -ssrc)
% We TESTED already in onebeam that abs(phi) < PHIBMAX

 global PHIBMAX;  % max phib (angle from equator) of our kernel

   global TBLK1 TBLK12t TBLK24tt TBLK24t TBLK4
   global R12 R24 RTOF;
   global RG12 RG24;

 global Ny1;
 global Ny2;
 global Nphi;

 global syy1 ;  % the out-projected screen coords of (recon pts -ssrc)
 global syy2 ;  % Set up in vecmap.m  Passed to here as globals; NOTE this
                % has size of ff, xx, yy, zz
 global syy0 ;  % coord in beam direction
 global rrsq ;  % 3D radius squared ; added Sep24_19

    syy1 = abs(syy1); 
    syy2 = abs(syy2); 
    phi  = abs(phi); %Recall K(y1,y2,phi) is symmetric in all 3 variables
           % and the TBLKPHI is only for nonnegative variables

% Recall y2 = 2*u2./(1-.999*u2) in mkTBLKPHI gives stretched y2 coord
       % thus u = y/(2+.999y) is the index in [0,1)
    u1 = syy1./(2+.999*syy1); % indices in [0,1) ; a matrix of size(syy1)
    u1 = (Ny1-1)*u1;  % indices in [0,Ny1-1)
    u2 = syy2./(2+.999*syy2); % indices in [0,1) ; a matrix of size(syy1)
    u2 = (Ny2-1)*u2;  % indices in [0,Ny2-1)

% Recall phi = PHIBMAX*(1-(1-u33)^2 )) in mkTBLKPHI gives stretched phi coord
    u33  = 1- sqrt(1-(phi/PHIBMAX)); % index in [0,1),  scalar
%   u33  = 1- sqrt(1-(phi/PHIBMAX)^2); % index in [0,1),  scalar ERROR july3_21
    u33  = (Nphi-1)*u33; % index in [0,Nphi-1), scalar

 f1 = round(u1) ; % integer part of u1, in [0,Ny1-1)
 f1 = 1+f1 ;                % Thus f1+1 is in [2,Ny1]

 f2 = round(u2) ; % integer part of u2, in [0,Ny2-1)
 f2 = 1+f2 ;

 f33 = floor(u33) ; % integer part of u33, in [0,Nphi-1) A SCALAR
 r33 = u33-f33;
 f33 = 1+f33 ;               % Thus f3+1 is in [2,Nphi]

    slice = (1-r33)*TBLK4(:,:,f33) + r33*TBLK4(:,:,1+f33);                ;

%  Jan23,07  I've figured out how to VECTORIZE this table lookup, by
%   realizing that the matrix TBLKPHI is stored as one long vector of length
%   Ny1*Ny2*Nph1;  the i,j,k matrix index has vector index 
%     i +Ny1*(j-1) +Ny1Ny2*(k-1) 
 
   Ny1Ny2 = Ny1*Ny2;

   u000 = f1 + Ny1*(f2-1) ;

        Hts = slice(u000);      % matrix of size(syy1),recon pts



%    cutoff for confidence wting.  NONE for f1 f2 f4

    
    % END of the vectorized interpolation

%   Hts gives the K values backprojected onto the matrix of recon pts; size(syy1) 
%   This Hts function was called by the main code "randPET3D"
   

