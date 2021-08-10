

% file randPET1july5_21B.m  (was randPET1july5_21.m  (was randPET124july3_21.m (was randPET124jan21_21.m (was randPET124jan4_21B.m  (was randPET124jan2_21B.m (was randPET124aug14_20B.m  (was randPET124aug13_20.m  (was randPET124may17_20.m  (was randPET124may5_20.m) (was ...mar31_20.m)  was ...nov8_19.m (randomized PET in 3D, with multigrid124)
  

%=================================================

   clear all
   close all   % close any old plots
    hold on    % hold on the new plots

%==================================================
%  We start the randomized PET reconstruction


   global PHIBMAX;   % max phi (angle from equator) of our kernel
   global xx yy zz;  % 3D matrices of x,y,z coords of the recon pts
   global syy1 syy2 syy0; % 3D matrices of out-projected screen coords 
                     % of the recon pts, for each onebeam
   global Ny1 Ny2 Nphi  ; % saved by mkTBLK.....
   global TBLK1  % saved by mkTBLK.....
%  ----------------------------------------

   global TOFu ;
   global N1 N2 N3 N4 N1234 Nsrc Nmax;
   global shift1 shift2 shift3 shift4;
   global count;
   global aa bb cc;

%  ----------------------------------------

   global xx yy zz ; % Feb26,07 I pass to vecmap these BIG input xx,yy,zz
                     % matrices of the xyz coords  of the recon pts

   t1 = clock;

   load svTBLK1 % pi/3 sin(phib) TBLK1 
 
   PHIBMAX    % print, from mkTBLK - should be pi/3
   TOFu = 10   % Print
%---------------------------

 N1=.5; N2=1.0; N3=.1; N4=1.0;
   N1234 = N1+N2+N3+N4;
   N1_N2_N3_N4_N1234  = [N1 N2 N3 N4 N1234]

%---------------------------

             %  Print the TBLKPHI parameters
    PHIBMAX_Ny1_Ny2_Nphi = [PHIBMAX,Ny1,Ny2,Nphi]

%----------------------------------------------------
 
%======================================================

        % the grid of reconstruction pts x,z,y
%     Nxx = 38;Nyy = 38; Nzz = 05;  % July5_21 
      Nxx = 38;Nyy = 38; Nzz = 00;  % July5_21 

  xx =   zeros(2*Nxx+1,2*Nyy+1,2*Nzz+1); 
  yy =   xx;
  zz =   xx;  % xx,yy,zz will contain the coords of this grid of recon pts
  syy0 = xx;  
  syy1 = xx;  
  syy2 = xx;  % the outprojected screen coords of these recon pts
  f1 =   xx;  %  CORRECTED version

%  Mar15,07  we fix the increments of the grid of recon pts
    dxx = 1.0;  dyy = 1.0;  dzz = 1.0;

  % set up the grid of x,y,z reconstruction pts in 3D

          for  kk = 1:2*Nzz+1
        for  jj = 1:2*Nyy+1
      for ii = 1:2*Nxx+1
  xx(ii,jj,kk) = (-Nxx+ii-1)*dxx; 
  yy(ii,jj,kk) = (-Nyy+jj-1)*dyy; 
  zz(ii,jj,kk) = (-Nzz+kk-1)*dzz; 
       end
         end
           end

  
  Nxx_Nyy_Nzz_dxx_dyy_dzz = [Nxx,Nyy,Nzz,dxx,dyy,dzz]
%==========================================================

  global aa bb cc Nsrc
  global count

%  aa = 1.0e-6; bb = 1.0e-6; cc = 1.0e-6; % half-dimensions of cube of sources
   aa = 08 ;bb = 08; cc = 02; % 


%  Nsrc = 1.0e2;
%  Nsrc = 1.0e3;
   Nsrc = 1.0e4;
%  Nsrc = 1.0e5;  % number of random sources to compute
%  Nsrc = 1.0e6;  % number of random sources to compute

   Nmax = N1234*Nsrc   % print
  aa_bb_cc_Nsrc =  [aa,bb,cc,Nsrc]   % print
  
   shift1 = [ 20 -14  00]; % for onebeam4cube %
   shift2 = [ 12  12  00]; % for onebeam4cube %
   shift3 = [ -14 20  00]; % for onebeam4cube %
   shift4 = [ 00  00  -10]; % for onebeam4cube %
   shift1_shift2_shift3_shift4 = [shift1 shift2 shift3 shift4] % print

   t1 = clock;
  

%========================================================

  count = 0;
 
              for ksrc = 1:Nmax 
  count = count + 1;

%  onebeam;  %This computes a SINGLE random source src inside the
     %  [-aa,aa]X[-bb,bb]X[-cc,cc] box  and returns ssrc,dir
     % where ssrc = src+TOFu*rand*dir is a pt shifted along the direction of the beam
%  onebeam1cube;
   onebeam4cubeMar20;



 % Now for each beam we BACK-PROJECT the kernel K onto the grid of RECON PTS
 % xx,yy,zz. This involves OUT-projecting these 3D pts onto the theta,phi 
 % detector screen to find their 2D parallel-beam screen coords syy1,syy.
 % (each of ff xx.....syy2 is vectorized as an 3D matrix) 
   
     
%  We already know that the phi angle of this single beam is  < PHIBMAX.

   phi = abs(asin(dir(3))); % make phi NONNEGATIVE since kernel K is symmetric
  
   vecmapsep24(ssrc,dir); %  Computes OUTPROJECTED screen coords syy0 syy1,syy2 of [(x,y,z)-ssrc]

          % and passes syy1 syy2 back as globals;
  
        % The syy1,syy2,phi have compressed values u1,u2,u3 in the
        % [0,1]X[0,1]X[0,1] cube. One then uses linear
        % interpolation on u3 and rounding on u1 and u2 in the TBLKPHI table to get
        %        K(syy1,syy2,phi)
        %  lookupK1(phi); % gives the BACKPROJECTION of this kernel.
        %  Recall that f starts out as zero

   f1 =   f1 + lookup1july5(phi); %    uses K1

           endfor    % end of ksrc loop over the Nsrc emission sources

   integral = dxx*dyy*dzz*sum(sum(sum(f1)));
   avintegral = integral/(1.6*Nsrc)  % print The 1.6 is needed july18_21

   mmmaxf1 = max(max(max(f1)))    % print
   mmminf1 = min(min(min(f1)))   % print July18_21 better
   minOmaxf1 = mmminf1/mmmaxf1   % print

   cpurecon = etime(clock,t1) % print

   
%%%%===================================================

    save svff.mat f1 PHIBMAX Ny1 Ny2 Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc Nmax aa bb cc TOFu N1 N2 N3 N4 N1234 mmmaxf1 mmminf1 minOmaxf1  cpurecon avintegral  % output string saved in svff.mat
          clear all
          close all
          load svff.mat

  mmmaxf1 = max(max(max(f1)))  % print
  ff1(:,:) = f1(:,:,Nzz+1);
  ff1 = (1/mmmaxf1)*ff1;


%%%%====================================================
     % Now some graphical output to visualize the reconstructed ff1 etc


%-----------------------
    hold on

  figure (2); mesh(ff1);
  colormap(jet); legend('ff1');
    hold on

  figure (3); plot(ff1);
  colormap(jet); legend('ff1');
    hold on

  figure (4); plot(min(0.2,ff1));
  colormap(jet); legend('ff1');
    hold on
 
  figure (5); contour(ff1);
  colormap(jet); colorbar; legend('ff1');
    hold on

  figure (6); image(60*ff1);
  colormap(jet); colorbar; 
    hold on

%    END of the reconstruction of the image ff on 4 cubes of recon pts

