
   %  file mkTBLK1july5_21.m  (was mkTBLK124july3_21.m (was mkTBLK124jan21_21.m  (was mkTBLK124jan17_21.m (was mkTBLK124jan6_21.m  (was mkTBLK124oct2_20.m  (was mkTBLK124aug14_20.m  (was mkTBLK124may5_20.m   ( was mkTBLK124mar31.m)

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

  PHIBMAX = pi/3

   Ny1 = 025; % test values to get SMALL svTBLK1
   Ny2 = 025; 
   Nphi =011;


%  Ny1 = 100; % values used in the paper
%  Ny2 = 100; 
%  Nphi =040;


   PHIBMAX_Ny1_Ny2_Nphi_NHu =  [PHIBMAX,Ny1,Ny2,Nphi,NHu] % print parameters

%---------------------------------------------


% Compressed grid points for our 3D matrix of the K(y1,y2,phi) kernel function
% We map [2000]X[2000]X[0,PHIBMAX] into the unit cube [0,1]X[0,1]X[0,1]

    TBLK1 = zeros(Ny1,Ny2,Nphi);
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

%---------------------------------------------
   
               for k3 = 1:Nphi 
   phi = pphi(k3); % in [0,PHIBMAX]; highly compressed for phi
     % near PHIBMAX, since our K has a sqrt(PHIBMAX-phi) type behavior there, 
     % and we want our table to allow accurate linear interpolation.
%====================================================     
%  NOTE: under table lookup in the routine lookupKphi.m, a point with    
%  screen coordinates sy1,sy2,phi will have indices in the [0,1) cube of 
%   u1 = sy1/(2+.999*sy1)
%   u2 = sy2/(2+.999*sy2)
%%%   u3 = 1- sqrt(1-(phi/PHIBMAX)^2)  %  ERROR July3_21
%   u3 = 1- sqrt(1-(phi/PHIBMAX))      %  CORRECTED July3_21

%=====================================================
    % WE BEGIN computing TBLK124 on the yy1,yy2 screen for each FIXED phi
       KK1 = zeros(Ny1,Ny2);

     % Feb9_2017 I do the phib = phi + vb^2 substitution for the integral
        % from phib = phi  to  phib = PHIBMAX, in order to handle the
        % 1/sqrt(PHIBMAX -phi) type singularity at phib near PHIBMAX. 

   dvb = sqrt(PHIBMAX+.00000001 - phi)/999; % the 999 is arbitrary, but 
      % is so large because of the large y1 y2 rays
      % The .00000001 is to avoid phib=phi and dividing by zero in
      %   vecHH1(phi,phib) below

               for kphib = 1:999
   vb = (kphib-.5)*dvb;  % I use Midpt rule for the integral, in order
                      % to avoid perhaps dividing by zero at the endpoints
   phib = phi + vb^2;

     
  etaphib = sin(phib);
  Kphib1 = etaphib*vecHH1july5(phi,phib); % blur radius 1.0
  % Returned is screen matrix of these K for fixed phi and phib

       KK1 = KK1 + 2*vb*dvb* Kphib1;
%================================================================
            end   % end of the kphib loop

   TBLK1(:,:,k3) = KK1 ;

                  end   % end of the k3 loop
 
   

%   Now the TBLK1(Ny1,Ny2,Nphi) is complete. To be used LATER for table
%     lookup to compute K(|syy1-y1*|,|syy2-y2*|,|phi|), where y1*,y*2 and syy1,syy2 are the
%     OUT-projected screen coords of the source and of the grid of xyz reconstr pts

    cpuTBLK124 = etime(clock,t0)

%%-------------------------------------------------------------
%    Now save the TBLK1 and various variables;

   save svTBLK1   PHIBMAX TBLHu TBLK1 Ny1 Ny2 Nphi NHu TBLHu 
    
%               We  now do output  visualization
   close all
   clear all
   load svTBLK1;
   
   PHIBMAX = PHIBMAX  % print
   
   mmmax1 = max(max(max(TBLK1)))   % print
   inv1 =        1/mmmax1;
   TBLK1 =       inv1*TBLK1;
   
   
         hold on
   TBLK1_1(:,:) = TBLK1(:,:,1);
   figure(2);  mesh(TBLK1_1); 
    legend('TBLK1_1')
      
         hold on
   figure(3);  plot(TBLK1_1) ; 
    legend('TBLK1_1')
      
         hold on
   TBLK1_9(:,:) = TBLK1(:,:,9);
   figure(4);  mesh(TBLK1_9); 
    legend('TBLK1_9')
      
         hold on
   figure(5);  plot(TBLK1_9) ; 
    legend('TBLK1_9')
      
         hold on
   TBLK1o1o(:,:) = TBLK1(:,1,:);  % plot does not work with 3 indices
   figure(6); plot(TBLK1o1o);
    legend('TBLK1o1o')
      
         hold on
    TBLK1o1op = TBLK1o1o' ; % plot does not work with " ' "
   TBLK1o1o(:,:) = TBLK1(:,1,:);
   figure(7); plot(TBLK1o1o');
    legend('TBLK1o1op')
      
         hold on
    figure(9);     plot(TBLHu);
    legend('TBLHu')

          hold off

%=================================================
%=================================================

