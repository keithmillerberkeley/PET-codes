

% file randPET124july24_21.m  (was randPET124feb2_21B.m  (was randPET124jan7_21B.m  (was randPET124jan5_21B.m (was randPET124jan4_21B.m  (was randPET124jan2_21B.m (was randPET124aug14_20B.m  (was randPET124aug13_20.m  (was randPET124may17_20.m  (was randPET124may5_20.m) (was ...mar31_20.m)  was ...nov8_19.m (randomized PET in 3D, with multigrid124)
  

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
   global TBLK1 TBLK2 TBLK4 TBLK1t TBLK2t TBLK12t TBLK24tt  % saved by mkTBLK.....
   global R12 R24 R1t fac;  % saved by mkTBLK.....
   global TOFu RG12 RG24 RG1t RG2t; 

   global xx yy zz ; % Feb26,07 I pass to vecmap these BIG input xx,yy,zz
                     % matrices of the xyz coords  of the recon pts

   global  syy0 syy1 syy2; % output matrices from vecmap of out-projected 
                     % screen coords syy1,syy2 of the recon pts
                     % for each single beam
   t1 = clock;


    load svTBLK124 % pi/6 sin(phib) TBLK1 TBLK2 TBLK4 etc R1t=20 R12=8 R24=16 fac=.8200,TOFu=10


             %  Print the TBLKPHI parameters
    PHIBMAX_Ny1_Ny2_Nphi = [PHIBMAX,Ny1,Ny2,Nphi]

   R1t = R1t  % from svTBLK124
   R12 = R12  % from svTBLK124
   R24 = R24  % from svTBLK124
%----------------------------------------------------
 
   TOFu = 10
%  TOFu = 03
   RG12 = 2.0*TOFu + 1.25*R12
   RG24 = 2.0*TOFu + 1.25*R24
   RG1t = 2.0*TOFu + 1.25*R1t

   fac= fac       %  May17_20

%======================================================

        % the grid of reconstruction pts x,z,y
      Nxx = 99;Nyy = 03; Nzz = 00;  %  Jan5_21

  xx =   zeros(2*Nxx+1,2*Nyy+1,2*Nzz+1); 
  yy =   xx;
  zz =   xx;  % xx,yy,zz will contain the coords of this grid of recon pts
  syy0 = xx;  
  syy1 = xx;  
  syy2 = xx;  % the outprojected screen coords of these recon pts
  f12t =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts
  f24tt =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts
  f4 =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts
  f1 =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts
  f1t =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts

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

%  aa = 1.0e-6; bb = 1.0e-6; cc = 1.0e-6; % dimensions of our box of sources
   aa = 4; bb = 4; cc = 16; % dimensions of our box of sources


%  Nsrc = 1.0e3;
%  Nsrc = 1.0e4;
   Nsrc = 1.0e5;  % number of random sources to compute
%%%   Nsrc = 1.0e6;  % number of random sources to compute

  aa_bb_cc_Nsrc =  [aa,bb,cc,Nsrc]   % print

%========================================================

  count = 0;
 
              for ksrc = 1:Nsrc 
  count = count + 1;

%  onebeam;  %This computes a SINGLE random source src inside the
     %  [-aa,aa]X[-bb,bb]X[-cc,cc] box  and returns ssrc,dir
     % where ssrc = src+3*dir is a pt shifted along the direction of the beam
   onebeam1cube;


 % Now for each beam we BACK-PROJECT the kernel K onto the grid of RECON PTS
 % xx,yy,zz. This involves OUT-projecting these 3D pts onto the theta,phi 
 % detector screen to find their 2D parallel-beam screen coords syy1,syy.
 % (each of ff xx.....syy2 is vectorized as an NxxXNyyXNphi matrix) 
   
     
%  We already know that the phi angle of this single beam is  < PHIBMAX.

   phi = abs(asin(dir(3))); % make phi NONNEGATIVE since kernel K is symmetric
  
   vecmapsep24(ssrc,dir); %  Computes OUTPROJECTED screen coords syy0 syy1,syy2 of [(x,y,z)-ssrc]

          % and passes syy1 syy2 back as globals;
  
        % The syy1,syy2,phi have compressed values u1,u2,u3 in the
        % [0,1]X[0,1]X[0,1] cube. One then uses 
        % trilinear interpolation on the TBLKPHI table to get
        %        K(syy1,syy2,phi)
        %  lookupK(phi); % gives the BACKPROJECTION of this kernel.
        %  Recall that f starts out as zero

   f1 =   f1 + lookup1aug14(phi); % uses  K12t truncated  
   f1t =   f1t + lookup1taug14(phi); % uses  K12t truncated  
   f12t =   f12t + lookup12taug14(phi); % uses  K12t truncated  
   f24tt=   f24tt + lookup24ttaug14(phi); % uses K24tt truncated
   f4 =   f4 + lookup4aug14(phi); % uses K4 NOT truncated 

           endfor    % end of ksrc loop over the Nsrc emission sources

   f124 = f12t+ fac*f24tt+ fac*fac*f4; % multigrid approx reconstruction of f1

   maxf124 = 1.0001*max(max(max(f124)))    % print
   imax124 = 1/maxf124;

   f1 = imax124*f1;
   f1t = imax124*f1t;
   f124 = imax124*f124;
   f12t = imax124*f12t;
   f24tt = imax124*f24tt;
   f4 =imax124*f4;

   cpurecon = etime(clock,t1) % print

   
%%%%===================================================

    save svff.mat f1 f1t f124 f12t f24tt f4 PHIBMAX Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc cpurecon  aa bb cc Nsrc  R1t R12 R24 RG12 RG24 TOFu fac ;  % output string saved in svff.mat
          clear all
          close all
          load svff.mat

% ff12t(:,:) = f12t(:,Nyy+1,:);
% ff24tt(:,:) = f24tt(:,Nyy+1,:);
%  ff4(:,:) = f4(:,Nyy+1,:);
% ff124(:,:) = f124(:,Nyy+1,:);
  
  ff1(:,:) = f1(:,:,Nzz+1);
  ff1t(:,:) = f1t(:,:,Nzz+1);
  ff124(:,:) = f124(:,:,Nzz+1);
  ff12t(:,:) = f12t(:,:,Nzz+1);
  ff24tt(:,:) = f24tt(:,:,Nzz+1);
  ff4(:,:) = f4(:,:,Nzz+1);


  ff24ttfill = fill2D2to1(ff24tt);

  ff4fill = fill2D4to1(ff4);

  ff124fill = ff12t + fac*ff24ttfill + fac*fac*ff4fill;


%%%%====================================================
     % Now some graphical output to visualize the reconstructed ff1 etc


%-----------------------
    hold on

  figure (2); plot(ff1);
  colormap(jet); legend('f1');
    hold on

  figure (3); plot(ff1t);
  colormap(jet); legend('f1t');
    hold on

  figure (4); plot(ff12t);
  colormap(jet); legend('f12t');
    hold on

  figure (5); plot(ff24tt);
  colormap(jet); legend('f24tt');
    hold on

  figure (6); plot(ff4);
  colormap(jet); legend('f4');
    hold on

  figure (7); plot(ff124);
  colormap(jet); legend('f124');
    hold on

  figure (8); plot(ff124fill);
  colormap(jet); legend('f124fill');
    hold on

  figure (9); plot(min(0.1,ff124fill));
  colormap(jet); legend('f124fill');
    hold on

  figure (10); plot(min(0.1,ff124));
  colormap(jet); legend('f124');
    hold on


%    END of the reconstruction of the image ff on a box of recon pts

