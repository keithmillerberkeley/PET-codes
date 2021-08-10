
%   file mkTBLK124jan29_21.m  (was mkTBLK124may5_20.m   ( was mkTBLK124mar31.m)

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

%    We begin the TBLKPHI computation

   t0 = clock;

   global yy1;      % Ny1XNy2
   global yy2;      % Ny1XNy2
   global PHIBMAX;
   global R1t R12 R24;

   PHIBMAX = 'pi/6'
   PHIBMAX = pi/6;

%  Ny1 = 100; % values used for Fig9 
%  Ny2 = 100; 
%  Nphi =040;

   Ny1 = 025;  % quick test values to produce a SMALL svTBLK124 file
   Ny2 = 025; 
   Nphi =011;

%  Ny1 = 040;  % these values probably OK for Fig9
%  Ny2 = 040; 
%  Nphi =010;
   
   R1t =20
   R12=08
   R24 =16  % always = 2*R12
%  R12=05
%  R24 =10
   fac = .8200; % tuned value for R12=8 RTOF=10 and PHIBMAX=pi/6
%  fac = .7141; % tuned value for R12=5 RTOF=10 and PHIBMAX=pi/6

   PHIBMAX_Ny1_Ny2_Nphi_NHu =  [PHIBMAX,Ny1,Ny2,Nphi,NHu] % print parameters

%---------------------------------------------


% Compressed grid points for our 3D matrix of the K(y1,y2,phi) kernel function
% We map [2000]X[2000]X[0,PHIBMAX] into the unit cube [0,1]X[0,1]X[0,1]

    TBLK1 = zeros(Ny1,Ny2,Nphi);
    TBLK2 = zeros(Ny1,Ny2,Nphi);
    TBLK4 = zeros(Ny1,Ny2,Nphi);
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

   rrr = sqrt(yy1.^2 + yy2.^2); % 2Dradius; to be used to cutoff at rrr > R12 or R24

%---------------------------------------------
   
               for k3 = 1:Nphi 
   phi = pphi(k3); % in [0,PHIBMAX]; highly compressed for phi
     % near PHIBMAX, since our K has a sqrt(PHIBMAX-phi) type behavior there, 
     % and we want our table to allow accurate linear interpolation.
     
%  NOTE: under table lookup in the routine lookupKphi.m, a point with    
%  screen coordinates sy1,sy2,phi will have indices in the [0,1) cube of 
%   u1 = sy1/(2+.999*sy1)
%   u2 = sy2/(2+.999*sy2)
%   u3 = 1- sqrt(1-(phi/PHIBMAX)^2)

%=====================================================
    % WE BEGIN computing TBLK124 on the yy1,yy2 screen for each FIXED phi
       KK1 = zeros(Ny1,Ny2);
       KK2 = zeros(Ny1,Ny2);
       KK4 = zeros(Ny1,Ny2);

       KK1t = zeros(Ny1,Ny2);
       KK2t = zeros(Ny1,Ny2);
       KK2tt = zeros(Ny1,Ny2);
       KK4tt = zeros(Ny1,Ny2);

     % Feb9_2017 I do the phib = phi + vb^2 substitution for the integral
        % from phib = phi  to  phib = PHIBMAX, in order to handle the
        % 1/sqrt(PHIBMAX -phi) type singularity at phib near PHIBMAX. 

   dvb = sqrt(PHIBMAX+.00000001 - phi)/999; % the 999 is large and arbitrary
      % The .00000001 is to avoid phib=phi and dividing by zero in
      %   vecHH(phi,phib) below

               for kphib = 1:999
   vb = (kphib-.5)*dvb;  % I use Midpt rule for the integral, in order
                      % to avoid perhaps dividing by zero at the endpoints
   phib = phi + vb^2;

     
%  pp = (PHIBMAX-phib)/PHIBMAX;  % this goes from 1 to 0
%  pp = sqrt(.0000001+pp);  % this goes from 1 to 0
%  etaphib = sin(phib)*pp;      %   Dec20_18
 
  etaphib = sin(phib); % value used for all Figs5-12

  Kphib1 = etaphib*vecHH1aug20(phi,phib); % blur radius 1.0
  Kphib2 = .125*etaphib*vecHH2aug20(phi,phib); % blur radius 2.0
  Kphib4 = .015625*etaphib*vecHH4aug20(phi,phib); % blur radius 4.0
  % Returned is screen matrix of these K for fixed phi and phib

       KK1 = KK1 + 2*vb*dvb* Kphib1;
       KK2 = KK2 + 2*vb*dvb* Kphib2;
       KK4 = KK4 + 2*vb*dvb* Kphib4;
%================================================================
            end   % end of the kphib loop

   KK1t = 0.5*(1-sign(-R12 +rrr )).*KK1;  % we cutoff at radius R12
   KK2t = 0.5*(1-sign(-R12 +rrr )).*KK2;  % we cutoff at radius R12
   KK2tt = 0.5*(1-sign(-R24 +rrr )).*KK2;  % we cutoff at radius R24=2*R12
   KK4tt = 0.5*(1-sign(-R24 +rrr )).*KK4;  % we cutoff at radius R24


   TBLK1(:,:,k3) = KK1 ;
   TBLK2(:,:,k3) = KK2 ;
   TBLK4(:,:,k3) = KK4 ;

   TBLK1t(:,:,k3) = KK1t ;
   TBLK2t(:,:,k3) = KK2t ;
   TBLK2tt(:,:,k3) = KK2tt ;
   TBLK4tt(:,:,k3) = KK4tt ;

                  end   % end of the k3 loop
 
   R12 = R12    %  print this is R12=8

   fac = fac   % Jan15_21 The tuned value is .8200 for R12=5 RTOF = 10


   TBLK12t = TBLK1t-fac*TBLK2t;  %
   TBLK24tt = fac*TBLK2tt-fac*fac*TBLK4tt;
   TBLK4    = fac*fac*TBLK4;


%   Now the TBLK124(Ny1,Ny2,Nphi) is complete. To be used LATER for table
%     lookup to compute K(syy1,syy2,phi), where syy1,syy2 are the "LOR-centered"
%     single-beam OUT-projected screen coords of the box of xyz reconstr pts

    cpuTBLK124 = etime(clock,t0)

%%-------------------------------------------------------------
%    Now save the TBLK and various variables;

   save svTBLK124   PHIBMAX TBLHu R12 R24 TBLK1 TBLK2 TBLK4 TBLK12t TBLK24tt Ny1 Ny2 Nphi NHu TBLHu  fac
    
%               We  now do output  visualization
   close all
   clear all
   load svTBLK124;
   
   mmmax1 = max(max(max(TBLK1))) % print
   mmmax2 = max(max(max(TBLK2)))
   mmmax4 = max(max(max(TBLK4)))
   
   inv1 = 1/mmmax1;
   TBLK1 = inv1*TBLK1;
   TBLK2 = inv1*TBLK2;
   TBLK4 = inv1*TBLK4;
   TBLK12t = inv1*TBLK12t;
   TBLK24tt = inv1*TBLK24tt;
   
   
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


         hold on

    figure(9);     plot(TBLHu);
    legend('TBLHu')
      hold on

          hold off

%=================================================
%=================================================

