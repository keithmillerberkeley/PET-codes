

   % file randPET1tjuly6_21.m  (was randPET1july5_21B.m  (was randPET1july5_21.m  (was randPET124july3_21.m (was randPET124jan21_21.m (was randPET124jan4_21B.m  (was randPET124jan2_21B.m (was randPET124aug14_20B.m  (was randPET124aug13_20.m  (was randPET124may17_20.m  (was randPET124may5_20.m) (was ...mar31_20.m)  was ...nov8_19.m (randomized PET in 3D, with multigrid124)
  

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
   global Ny1 Ny2 Nphi  ; % saved by mkTBLK1t.....
   global TBLK1t  % saved by mkTBLK1t.....
%  ----------------------------------------

   global R1t % saved by mkTBLK1t.....
   global TOFu LG1t;
   global N1 N2 N3 N4 N1234 Nsrc Nmax;
   global shift1 shift2 shift3 shift4;
   global count;
   global aa bb cc;

%  ----------------------------------------

   global xx yy zz ; % Feb26,07 I pass to vecmap these BIG input xx,yy,zz
                     % matrices of the xyz coords  of the recon pts

   t1 = clock;

   load svTBLK1t % pi/6 sin(phib) TBLK1t 

  

%  R1t = R1t    % =10 from MKTBLK
%  TOFu = 20   % Print TEMPOR 

%  R1t = R1t    % =20 from MKTBLK July7_21
%  TOFu = 10   % Print july19_21

%  R1t = R1t    % =10 from MKTBLK Aug3_21
%  TOFu = 20   % 

   R1t = R1t    % =10 from MKTBLK Aug3_21
   TOFu = 30   % Print % use with  LG1t = 0.0R1t +1.0TOFu

%---------------------------
   LG1t = 0.0*R1t + 1.0*TOFu  % TEMPOR aug3_21 to test the "2.0"
%  LG1t = 1.25*R1t + 0.0*TOFu  % TEMPOR aug3_21 to test the "2.0"
%  LG1t = 1.25*R1t + 2.0*TOFu

%---------------------------

   N1=.5; N2=1.0; N3=.1; N4=1.0;
   N1234 = N1+N2+N3+N4;
   N1_N2_N3_N4_N1234  = [N1 N2 N3 N4 N1234]

%---------------------------

             %  Print the TBLK1t parameters
    PHIBMAX_Ny1_Ny2_Nphi = [PHIBMAX,Ny1,Ny2,Nphi]

%----------------------------------------------------
 
%======================================================

        % the grid of reconstruction pts x,z,y
      Nxx = 38;Nyy = 38; Nzz = 00;  % July6_21 and July19_21
%     Nxx = 38;Nyy = 38; Nzz = 05;  % July6_21  Aug2_21

  xx =   zeros(2*Nxx+1,2*Nyy+1,2*Nzz+1); 
  yy =   xx;
  zz =   xx;  % xx,yy,zz will contain the coords of this grid of recon pts
  syy0 = xx;  
  syy1 = xx;  
  syy2 = xx;  % the outprojected screen coords of these recon pts
  f1t =   xx;  %  CORRECTED version

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

   aa = 08 ;bb = 08; cc = 02; % 


%   Nsrc = 1.0e2;
%   Nsrc = 1.0e3;
    Nsrc = 1.0e4;
%   Nsrc = 1.0e5;  % number of random sources to compute
%%%%%    Nsrc = 1.0e6;  % number of random sources to compute

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
        % interpolation on u3 and rounding on u1 and u2 in the TBLK1t table to get
        %        K1t(syy1,syy2,phi)
        %  lookupK1t(phi); % gives the BACKPROJECTION of this kernel.
        %  Recall that f starts out as zero

   f1t =   f1t + lookup1tjuly6(phi); %    uses K1t

           endfor    % end of ksrc loop over the Nsrc emission sources

   integral1t = dxx*dyy*dzz*sum(sum(sum(f1t)));
   avintegral1t = integral1t/Nsrc  % print

   mmmaxf1t = max(max(max(f1t)))    % print
   mmminf1t = min(min(min(f1t)))
   minOmaxf1t = mmminf1t/mmmaxf1t

   cpurecon = etime(clock,t1) % print

   
%%%%===================================================

    save svff.mat f1t PHIBMAX R1t LG1t Ny1 Ny2 Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc Nmax aa bb cc TOFu N1 N2 N3 N4 N1234 mmmaxf1t mmminf1t minOmaxf1t cpurecon avintegral1t  % output string saved in svff.mat
          clear all
          close all
          load svff.mat

  ff1t(:,:) = f1t(:,:,Nzz+1);
  ff1t = (1/mmmaxf1t)*ff1t;
   intneg = sum(sum(sum(min(0.0,f1t))))
   int    = sum(sum(sum(f1t)))
   intnegOint = intneg/int



%%%%====================================================
     % Now some graphical output to visualize the reconstructed ff1 etc


%-----------------------
    hold on

  figure (2); mesh(ff1t);
  colormap(jet); legend('ff1t');
    hold on

  figure (3); plot(ff1t);
  colormap(jet); legend('ff1t');
    hold on

  figure (4); plot(min(0.2,ff1t));
  colormap(jet); legend('ff1t');
    hold on
 
  figure (5); contour(ff1t);
  colormap(jet); colorbar; legend('ff1t');
    hold on

  figure (6); image(60*ff1t);
  colormap(jet); colorbar; 
    hold on

  figure (7); plot(min(.2,ff1t));
  colormap(jet); legend('ff1t');
    hold on

%    END of the reconstruction of the image ff on 4 cubes of recon pts

