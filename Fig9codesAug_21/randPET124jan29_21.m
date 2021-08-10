

%     file randPET124jan29_21.m  (was randPET124may5_20.m (was ...mar31_20.m)  was ...nov8_19.m (randomized PET in 3D, with multigrid124)
  

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
   global Ny1 Ny2 Nphi  ; % saved by mkTBLK124.....
   global TBLK1 TBLK1t TBLK12t TBLK24tt TBLK4 TBLK2  % saved by mkTBLK124.....
   global R12 R24 ;  % saved by mkTBLK124.....
   global RTOF RG12 RG24 ; 
   global N1 N2 N3 N4 N1234 Nsrc Nmax;
   global shift1 shift2 shift3 shift4;
   global aa bb cc; 
   global xx yy zz;
   global count
   global fac     % saved by mkTBLK124

       % NOTE 12t uses fac; 24tt uses fac and fac*fac; 4 uses fac*fac

 load svTBLK124  % pi/6 sin(phib) TBLK1 etc all 100X100X40 r12=5 and fac=.7141 tuned for RTOF = 10



             %  Print the TBLKPHI parameters
    PHIBMAX_Ny1_Ny2_Nphi = [PHIBMAX,Ny1,Ny2,Nphi]

   R12 = R12  % from svTBLKPHI
   R24 = R24  % from svTBLKPHI
%----------------------------------------------------
 
   RTOF = 10
   RG12 = 2.0*RTOF + 1.25*R12
   RG24 = 2.0*RTOF + 1.25*R24
   fac= fac       %  May17_20

   N1=.5; N2=1.0; N3=.1; N4=1.0;
   N1234 = N1+N2+N3+N4;
   N1_N2_N3_N4_N1234  = [N1 N2 N3 N4 N1234]

%  Nsrc = 1.0e2;
%  Nsrc = 1.0e3;
   Nsrc = 1.0e4;
%  Nsrc = 1.0e5;  % number of random sources to compute
%%%%%  Nsrc = 1.0e6;  % number of random sources to compute.  (value for Fig9)
%  Nsrc = 1.0e7;  % number of random sources to compute

   Nmax = N1234*Nsrc

   shift1 = [ 20 -14  00]; % for onebeam4cube %
   shift2 = [ 12  12  00]; % for onebeam4cube %
   shift3 = [ -14 20  00]; % for onebeam4cube %
   shift4 = [ 00  00  -20]; % for onebeam4cube %
   shift1_shift2_shift3_shift4 = [shift1 shift2 shift3 shift4]

   t1 = clock;

%======================================================

        % the BOX of reconstruction pts x,z,y
%     Nxx = 02;Nyy = 02; Nzz = 00;
      Nxx = 38;Nyy = 38; Nzz = 00; % july23_21 Nxx should be EVEN for "fill4to1"
%     Nxx = 39;Nyy = 39; Nzz = 00;
%     Nxx = 48;Nyy = 48; Nzz = 10;

  xx =   zeros(2*Nxx+1,2*Nyy+1,2*Nzz+1); 
  yy =   xx;
  zz =   xx;

%%%%%%-------------------------------------
%  global xx yy zz ; % Feb26,07 I pass to vecmap these BIG input xx,yy,zz
                     % matrices of the xyz coords  of the recon pts
%  global  syy0 syy1 syy2; % output matrices from vecmap of out-projected 
                     % screen coords syy1,syy2 of the recon pts
                     % for each single beam
% global aa bb cc Nsrc
% global count
%%%%%%-------------------------------------

%  aa = .000001; bb = .000001; cc = .000001; % dimensions of our box of sources
%  aa = 05; bb = 05; cc = 02; % dimensions of our box of sources
   aa = 08; bb = 08; cc = 02; % dimensions of our box of sources

  aa_bb_cc_Nsrc_Nmax = [aa bb cc Nsrc Nmax]  % print

   t1 = clock;

%======================================================

  xx =   zeros(2*Nxx+1,2*Nyy+1,2*Nzz+1); 
  yy =   xx;
  zz =   xx;  % xx,yy,zz will contain the coords of this box of recon pts
  syy0 = xx;  
  syy1 = xx;  
  syy2 = xx;  % the outprojected screen coords of these recon pts
  f1 =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts
  f12t =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts
  f24tt =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts
  f4 =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts
  f2 =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts

%  Mar15,07  we fix the increments of the box of recon pts
    dxx = 1.0;  dyy = 1.0;  dzz = 1.0;

  % set up the BOX of x,y,z reconstruction pts in 3D

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

%  aa = .000001; bb = .000001; cc = .000001; % dimensions of our box of sources
%  aa = 05; bb = 05; cc = 02; % dimensions of our box of sources
   aa = 08; bb = 08; cc = 02; % dimensions of our box of sources

%========================================================

  count = 0;
 
              for ksrc = 1:Nmax 
  count = count + 1;

%  onebeam;  %This computes a SINGLE random source src inside the
     %  [-aa,aa]X[-bb,bb]X[-cc,cc] box  and returns ssrc,dir
     % where ssrc = src+3*dir is a pt shifted along the direction of the beam
%  onebeam1cube;
%  onebeam2cube;
   onebeam4cubeMar20


 % Now for each beam we BACK-PROJECT the kernel K onto the BOX of RECON PTS
 % xx,yy,zz. This involves OUT-projecting these 3D pts onto the theta,phi 
 % detector screen to find their 2D parallel-beam screen coords syy1,syy.
 % (each of ff xx.....syy2 is vectorized as an (2*Nxx+1)X(2*Nyy+1)X(2*Nzz+1)matrix) 
   
     
%  We already know that the phi angle of this single beam is  < PHIBMAX.

   phi = abs(asin(dir(3))); % make phi NONNEGATIVE since kernel K is symmetric
  
   vecmapsep24(ssrc,dir); %  Computes OUTPROJECTED screen coords syy0 syy1,syy2

          % and passes syy1 syy2 back as globals;
  
        % The syy1,syy2,phi have compressed values u1,u2,u3 in the
        % [0,1]X[0,1]X[0,1] cube. One then uses 
        % trilinear interpolation on the TBLKPHI table to get
        % K(syy1,syy2,phi)

        %  lookupKphi(phi); % gives the BACKPROJECTION of this kernel.
        %  Recall that ff starts out as zero

   f1 =   f1 + lookup1nov8(phi); % uses K1 NOT truncated
%  f1t =   f1t + lookup1tnov8(phi); % uses K1  truncated at r<R1t
   f12t =   f12t + lookup12tnov8(phi); % uses  K12t truncated at r<R12 
   f24tt=   f24tt + lookup24ttnov8(phi); % uses K24tt truncated at r<R24=2*R12
   f4 =   f4 + lookup4nov8(phi); % uses K4 NOT truncated 
   f2 =   f2 + lookup2nov8(phi); % uses K2 NOT truncated

           endfor    % end of ksrc loop over the Nsrc emission sources

   f124 = f12t+f24tt+f4; % multigrid approx reconstruction of f1
    cpurecon = etime(clock,t1) % print
   
%%%%===================================================

    save svff.mat f124 f12t  f24tt f1 f4 f2 PHIBMAX Ny1 Ny2 Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc cpurecon  aa bb cc Nsrc R12 R24 RG12 RG24 RTOF N1 N2 N3 N4 N1234 shift1 shift2 shift3 shift4 fac;  % output string saved in svff.mat

          clear all
          close all
          load svff.mat

% Now take  2D (x,0,z) slices, ff1 of f1 etc
  ff1(:,:) = f1(:,:,Nzz+1);
  ff124(:,:) = f124(:,:,Nzz+1);
  ff12t(:,:) = f12t(:,:,Nzz+1);
  ff24tt(:,:) = f24tt(:,:,Nzz+1);
   ff2(:,:) = f2(:,:,Nzz+1);
   ff4(:,:) = f4(:,:,Nzz+1);

  ff24ttfill = fill2D2to1(ff24tt);

  ff4fill = fill2D4to1(ff4);

  ff124fill = ff12t + ff24ttfill + ff4fill;


%%%%====================================================
     % Now some graphical output to visualize the reconstructed ff1 etc

  max1 = 1.01*max(max(ff1))
  inv1 = 1/max1;
  ff1 = inv1*ff1;
  ff124 = inv1*ff124;
  ff12t = inv1*ff12t;
  ff24tt = inv1*ff24tt;
  ff4 = inv1*ff4;
  ff2 = inv1*ff2;


  ff124fill = inv1*ff124fill;
  ff24ttfill = inv1*ff24ttfill;
  ff4fill = inv1*ff4fill;

  diff1m124fill = ff1-ff124fill;


%-----------------------
    hold on
  figure (2); mesh(ff1);
  colormap(jet); legend('ff1'); 
  % But Note: print -depsc2 will allow annotation of legend&plot, But NOT of legend&mesh. See view2Df124jan29_21.m and mk2X2newFig9jan31_21.m for that.
    hold on

  figure (3); mesh(ff124fill);
% figure (3); mesh(ff124);
% figure (3); mesh(min(.10,ff1));
% colormap(jet); legend('ff124');
  colormap(jet); legend('ff124fill');
    hold on

  figure (4); mesh(ff12t);
% figure (4); mesh(min(.10,ff1));
  colormap(jet); legend('ff12t');
    hold on

  figure (5); mesh(ff24tt);
% figure (5); mesh(min(.10,ff1));
  colormap(jet); legend('ff24tt');
    hold on

  figure (6); mesh(ff4fill);
% figure (6); mesh(ff4);
% figure (6); mesh(min(.10,ff1));
% colormap(jet); legend('ff4');
  colormap(jet); legend('ff4fill');
    hold on


     hold off

%    END of the reconstruction of the image ff on a box of recon pts

