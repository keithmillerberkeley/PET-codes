
 %  file mkTBLK124july14_21.m  (was mkTBLK124jan17_21.m (was mkTBLK124jan6_21.m  (was mkTBLK124oct2_20.m  (was mkTBLK124aug14_20.m  (was mkTBLK124may5_20.m   ( was mkTBLK124mar31.m)

%=================================================
   clear all
   close all
%==================================================
%  BEGIN computation by computing a table of NHu = 20001 values of Hu(u)
%    = H(t), where u = 20000*t*t/(4+t*t) on equispaced u values in [0,1].
%   We will linearly interpolate our H(u) values from this table.

 global TBLHu  % we make this table global (to be used by vecKphi.m)
 global NHu  % # of  TBLHu entries, usually 2001

 NHu = 2001;
 du = 1/(NHu-1);
 u = [0:du:1];
 TBLHu = zeros(1,NHu);

   for j = 1:NHu-1
 t = 4*u(j)/(1-u(j));
  TBLHu(j) = H(t)*(1 +t^2);  % 2D exp(.5*r^2);
    end
  TBLHu(NHu) = -1.0; % this is the value at t = infinity for H


%  end of table setup for Hu(u) to compute H(t) by lookup
%===================================================
%===================================================

%    We begin the TBLK computation  for f124fill

   t0 = clock;

   global yy1;      % Ny1XNy2
   global yy2;      % Ny1XNy2
   global PHIBMAX;
   global R12 R24 R1t fac;

   PHIBMAX = pi/6
   R1t = 20
%  R12=08
%  R24=16
%  fac = .9053  % results in fton = -.20132  with Nsrc=1.0e3 with TOFu=10
        % This is the zeroXsection fac= .9053 for R12=8
%  fac = .8500  % results in fton = -.05846  with Nsrc=1.0e3 with TOFu=10
%  fac = .8200 % results in fton = -.00511 OR +,00533  with Nsrc=1.0e3 with TOFu=10
        % Thus fac = .8200 is our TUNED value for R12=8, TOFu=10, pi/6

   R12 = 08   
   R24 = 16
   fac = .6120   % the zXs value % but this results in fton = -.002008 or 
        % -.00863 or+.00418 with Nsrc=1.0e3; OR = -.001244 with Nsrc=1.0e4
        % Thus fac = .6012 is our TUNED value for R12=8, TOFu=3, pi/6

%  R12 = 05   
%  R24 = 10
% fac = .82012   % the zXs value % results in fton = -.55975
% fac = .71410   % the ERROR value % results in fton = -.01349
% fac = .70500   % results in fton = +.02464
% fac = .71000   % results in fton = +.0115 or +.0200
% fac = .71100   % results in fton = -.00662 or = -.0128 with fton cutoff at r>6
                 % results in fton = -.00246 with fton cutoff at r>4
         % Thus fac = .7110 is our TUNED value for R12=5, TOFu=10, pi/6

%  R12 = 05
%  R24 = 10
%  fac = .82012  % the zXs value % results in fton = -.4384; with fton cutoff at r>4
%  fac = .6000   % results in fton = -.0666
%  fac = .5500   % results in fton = +.0202
%  fac = .5876   % results in fton = +.00655 (I used testsecant.m)
         % Thus fac = .5876 is our TUNED value for R12=5, TOFu=3, pi/6
   
%%%%%%%%%%%%%%%%     July16_21continue to find TUNED facs

%  PHIBMAX = atan(.8*15/40)   % with cbm and QZ="NOT=1.0" and Nsrc=1.0e3
%  R12 = 08
%  R24 = 16     %  for TOFu=10 and cbm
%  fac = .8209  %  from Fig11AB  % results in fton = +.005392 or +.01228
%  fac = .8220  %  fton = +.1993
%  fac = .8300  %  fton = -.0108
%  fac = .8257  %  fton = -.01965
%  fac = .8200  %  fton = +.00910
          % Thus fac = .8200 is ourTUNED value for R12=8, TOFu =10, "NOT=1.0"

%  PHIBMAX = atan(15/40) % with cbm and "QZ=1.0" and Nsrc=1.0e3
%  R12 = 08
%  R24 = 16     %  for TOFu=10 and cbm and Fig11C
%  fac = .8269  % from Fig11C, results fton= -.00880 with r>4 cutoff
%  fac = .8269  % from Fig11C, results fton= -.01218 with r>6 cutoff
%  fac = .8200  %              results fton= -.00978 with r>6 cutoff
%  fac = .8300  %              results fton= -.00166 with r>6 cutoff
          % Thus fac = .8300 is ourTUNED value for R12=8, TOFu =10, "YES=1.0"

%%%%% 

%  Ny1= 100
%  Ny2= 100
%  Nphi = 40

   Ny1= 25
   Ny2= 25
   Nphi = 11

   PHIBMAX_Ny1_Ny2_Nphi_NHu =  [PHIBMAX,Ny1,Ny2,Nphi,NHu] % print parameters

%---------------------------------------------


% Compressed grid points for our 3D matrix of the K(y1,y2,phi) kernel function
% We map [2000]X[2000]X[0,PHIBMAX] into the unit cube [0,1]X[0,1]X[0,1]

    TBLK1 = zeros(Ny1,Ny2,Nphi);
    TBLK2 = zeros(Ny1,Ny2,Nphi);
    TBLK4 = zeros(Ny1,Ny2,Nphi);
    TBLK1t = zeros(Ny1,Ny2,Nphi);
    TBLK12t = zeros(Ny1,Ny2,Nphi);
    TBLK24tt = zeros(Ny1,Ny2,Nphi);
    du1 = 1/(Ny1-1);
    du2 = 1/(Ny2-1);
    du3 = 1/(Nphi-1);
   u1=[0:du1:1];  
   u2=[0:du2:1];
   u3=[0:du3:1];

 y1 = 2*u1./(1-.999*u1); % vector of y1  values in our 3D matrix; in [0,2000]
 y2 = 2*u2./(1-.999*u2); % vector of y2  values in our 3D matrix; in [0,2000]
 pphi = PHIBMAX*(1-(1-u3).^2 ); % the phi values in our matrix; in [0,PHIBMAX]

              for k1 = 1:Ny1
           for k2 = 1:Ny2
   yy1(k1,k2) = y1(k1); 
   yy2(k1,k2) = y2(k2); 
           end
              end
% Now yy1,yy2 are Ny1XNy2 matrices of y1,y2 coords on 2D screen for
%    each fixed phi

   rrr = sqrt(yy1.^2 + yy2.^2); % 2Dradius; to be used to cutoff at rrr= R12 or R24

%---------------------------------------------
   
               for k3 = 1:Nphi 
   phi = pphi(k3); % in [0,PHIBMAX]; highly compressed for phi
     % near PHIBMAX, since our K has a sqrt(PHIBMAX-phi) type behavior there, 
     % and we want our table to allow accurate linear interpolation.
     
%  NOTE: under table lookup in the routine lookupKphi.m, a point with    
%  screen coordinates sy1,sy2,phi will have indices in the [0,1) cube of 
%   u1 = sy1/(2+.999*sy1)
%   u2 = sy2/(2+.999*sy2)
%   u3 = 1- sqrt(1-(phi/PHIBMAX))   % CORRECT
%   u3 = 1- sqrt(1-(phi/PHIBMAX)^2)  % ERROR found, July2021

%=====================================================
    % WE BEGIN computing TBLK124 on the yy1,yy2 screen for each FIXED phi
       KK1 = zeros(Ny1,Ny2);
       KK2 = zeros(Ny1,Ny2);
       KK4 = zeros(Ny1,Ny2);

       KK12t = zeros(Ny1,Ny2);
       KK24tt = zeros(Ny1,Ny2);

     % Feb9_2017 I do the phib = phi + vb^2 substitution for the integral
        % from phib = phi  to  phib = PHIBMAX, in order to handle the
        % 1/sqrt(PHIBMAX -phi) type singularity at phib near PHIBMAX. 

   dvb = sqrt(PHIBMAX+.00000001 - phi)/399; % the 999 is large and arbitrary
%  dvb = sqrt(PHIBMAX+.00000001 - phi)/999; % the 999 is large and arbitrary
      % The .00000001 is to avoid phib=phi and dividing by zero in
      %   vecHH(phi,phib) below

               for kphib = 1:399
%              for kphib = 1:999
   vb = (kphib-.5)*dvb;  % I use Midpt rule for the integral, in order
                      % to avoid perhaps dividing by zero at the endpoints
   phib = phi + vb^2;

     
  etaphib = sin(phib);
  Kphib1 = etaphib*vecHH1aug20(phi,phib); % blur radius 1.0
  Kphib2 = etaphib*vecHH2aug20(phi,phib); % blur radius 2.0
  Kphib4 = etaphib*vecHH4aug20(phi,phib); % blur radius 4.0
  % Returned is screen matrix of these K for fixed phi and phib

       KK1 = KK1 + 2*vb*dvb* Kphib1;
       KK2 = KK2 + 2*vb*dvb* Kphib2;
       KK4 = KK4 + 2*vb*dvb* Kphib4;
%================================================================
            end   % end of the kphib loop

   TBLK1(:,:,k3) = KK1 ;
   TBLK2(:,:,k3) = .125*KK2 ;
   TBLK4(:,:,k3) = .015625*KK4 ;

                  end   % end of the k3 loop
 
   TBLK1t = 0.5*(1-sign(-R1t -.000001 +rrr )).*TBLK1;                 % we cutoff at radius R1t
   TBLK12t = 0.5*(1-sign(-R12 -.000001 +rrr )).*(TBLK1 -fac*TBLK2);  % we cutoff at radius R12
   TBLK24tt = 0.5*(1-sign(-R24 -.000001 +rrr )).*(TBLK2 -fac*TBLK4);  % we cutoff at radius R24
   
%% NOTE later in randPET124;   f124 =  f12t + fac*f24tt + fac*fac*f4;

%   Now the TBLK124(Ny1,Ny2,Nphi) is complete. To be used LATER for table
%     lookup to compute K(|syy1-y1*|,|syy2-y2*|,|phi|), where y1*,y*2 and syy1,syy2 are the
%     OUT-projected screen coords of the source and of the grid of xyz reconstr pts

    cpuTBLK124 = etime(clock,t0)

%    STOP HERE 8:00pm aug14
%    STOP HERE 8:00pm aug14
%%-------------------------------------------------------------
%    Now save the TBLK124 and various variables;

   save svTBLK124   PHIBMAX TBLHu R1t R12 R24 TBLK1 TBLK2 TBLK4 TBLK1t TBLK12t TBLK24tt Ny1 Ny2 Nphi NHu TBLHu fac
    
%               We  now do output  visualization
   close all
   clear all
   load svTBLK124;
   
   mmmax1 = max(max(max(TBLK1)))   % print
   mmmax2 = max(max(max(TBLK2)))   % print
   mmmax4 = max(max(max(TBLK4)))   % print
   
   inv1 =        1/mmmax1;
   TBLK1 =       inv1*TBLK1;
   TBLK2 =       inv1*TBLK2;
   TBLK4 =       inv1*TBLK4;
   TBLK12t =     inv1*TBLK12t;
   TBLK24tt =    inv1*TBLK24tt;
   TBLK4 =       inv1*TBLK4;
   
   
         hold on
   figure(2);  mesh(TBLK1(:,:,1) ); 
%  figure(2);  plot(TBLK1(:,:,8) ); 
    legend('TBLK1')
      
         hold on
   figure(3);  mesh(TBLK2(:,:,1) ); 
%  figure(3);  plot(TBLK2(:,:,8) );
    legend('TBLK2')

         hold on
   figure(4);  mesh(TBLK4(:,:,1) ); 
%  figure(4);  plot(TBLK4(:,:,8) );
    legend('TBLK4')

         hold on
   figure(5);  mesh(TBLK12t(:,:,1) ); 
%  figure(5);  plot(TBLK12t(:,:,8) ); 
    legend('TLBK12t')

         hold on
   figure(6);  mesh(TBLK24tt(:,:,1) ); 
%  figure(6);  plot(TBLK24tt(:,:,8) ); 
    legend('TBLK24tt')

      
         hold on
    figure(9);     plot(TBLHu);
    legend('TBLHu')

          hold off

%=================================================
%=================================================

