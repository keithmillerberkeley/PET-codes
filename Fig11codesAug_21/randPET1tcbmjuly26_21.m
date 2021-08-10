
   %  file randPET1tcbmjuly26_21.m  (was randPET1tcbmfeb5_21.m  (was randPET124cbmjan24_21.m  (was randPET124cbmjan23_21.m  (wasrandPETf124cbmjan21_21.m  (was randPETf124cbmjan19_21.m (was randPETf1tcbmjan16_21.m  (was randPET1tcbmmar8_20 .m  (randomized PET in 3D)
  

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
   global TBLK1t  % saved by mkTBLK.....
   global R1t ; %saved by mkTBLK
   global TOFu;
   global xx yy zz ; % Feb26,07 I pass to vecmap these BIG input xx,yy,zz
               % matrices of the xyz coords  of the recon pts
   global  syy0 syy1 syy2; % output matrices from vecmap of out-projected 
               % screen coords syy1,syy2 of the recon pts for each single beam
   global aa bb cc Nsrc 
   global count
   global RG1t

   load svTBLK1t

             %  Print the TBLK parameters; from mkTBLK
    PHIBMAX_Ny1_Ny2_Nphi = [PHIBMAX,Ny1,Ny2,Nphi]
   R1t  % print
  
%===================================================

   TOFu = 20   % Jan16_21 and Jan22_21 andJan23_21
   RG1t = 2.0*TOFu + 1.25*R1t

   ttt1 = clock;

%  Nsrc = 1.0e3  ; 
   Nsrc = 1.0e4  ; 
%  Nsrc = 1.0e5  ; 
%  Nsrc = 1.0e6  ; 
%%%%   Nsrc = 1.0e7  ; 

%======================================================

   aa = 02; bb = 02; cc = 45; % dimensions of our box of sources 

  aa_bb_cc_Nsrc =  [aa,bb,cc,Nsrc]   % print


        % the grid of reconstruction pts x,z,y
     Nxx = 38;Nyy = 00; Nzz = 38;

  xx =   zeros(2*Nxx+1,2*Nyy+1,2*Nzz+1); 
  yy =   xx;
  zz =   xx;  % xx,yy,zz will contain the coords of this box of recon pts
  syy0 = xx;  
  syy1 = xx;  
  syy2 = xx;  % the outprojected screen coords of these recon pts
  f1t =   xx;  % zero the initial RECONSTRUCTED IMAGE on these recon pts

%  Mar15,07  we fix the increments of the box of recon pts
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

  Nxx_Nyy_Nzz_dxx_dyy_dzz = [Nxx,Nyy,Nzz,dxx,dyy,dzz]  % print
%==========================================================

   RR = 40                  % print  radius of scanner
   LL = 15  % print half-length of scanner
%  LOSW, desired halflength of scan window = RR for this simulation
   LOSW = RR  % Feb6,7  and Mar8
   LLL = (LOSW + LL +R1t*sec(PHIBMAX))   % print half-length center will travel in Nmax counts
   DZ  = 2*LLL/Nsrc     % z translation per count

     % Thus center of scanner is at Cz = -LLL+DZ*count; and scanner detects all z in [Cz-LL,Cz+LL]
     % and we want z1,z2 inside. That is, max([abs(z1-Cz),abs(z2-Cz)])<2*LL
%  ERROR discovered Feb7:  % we want z1,z2 inside. That is, max([abs(z1-Cz),abs(z2-Cz)])<LL

     count = 0;
  rejecLL = 0;
  rejecD  = 0;

              for ksrc = 1:Nsrc

%  onebeam;  %This computes a SINGLE random source src inside the
     %  [-aa,aa]X[-bb,bb]X[-cc,cc] box  and returns ssrc,dir
     % where ssrc = src+3*dir is a pt shifted along the direction of the beam
   onebeam1cubeMar8_20;  %  rotation by 45degrees 

 % Now for each beam we BACK-PROJECT the kernel K onto the BOX of RECON PTS
 % xx,yy,zz. This involves OUT-projecting these 3D pts onto the theta,phi 
 % detector screen to find their 2D parallel-beam screen coords syy1,syy.
 % (each of ff xx.....syy2 is vectorized as an 2*Nxx+1,2*Nyy+1,2*Nzz+1 matrix) 
   
     
%  We already know that the phi angle of this single beam is  < PHIBMAX.

   phi = abs(asin(dir(3))); % make phi NONNEGATIVE since kernel K is symmetric
  
   vecmap(ssrc,dir); %  Computes OUTPROJECTED screen coords syy0 syy1,syy2
          % of this beam from the box of recon pts xx,yy,zz 
          % and passes syy1 syy2 back as globals;
  
        % The syy1,syy2,phi have compressed values u1,u2,u3 in the
        % [0,1]X[0,1]X[0,1] cube. One then uses 
        %  TBLKPHI table to get
        % K(syy1,syy2,phi)

        %  lookupKphi(phi); % gives the BACKPROJECTION of this kernel.
        %  Recall that ff starts out as zero
%%=================================================

       sx=ssrc(1); sy=ssrc(2); sz=ssrc(3);
       dx=dir(1); dy=dir(2); dz=dir(3);

       C =  (-RR^2+sx^2+sy^2);
       B = sx*dx+sy*dy;
       A = dx^2+dy^2;

       D = (B^2-A*C);         % print
              if (D<0)
       rejecD = rejecD+1;

              else
       D = sqrt(D);         % print
       t1 = (-B-D)/A;
       t2 = (-B+D)/A;
       z1 = dz*t1 + sz;    % print
       z2 = dz*t2 + sz;   % print
       z2mz1 = z2-z1;
                    QZ = 2*LL/(2*LL-abs(z2-z1)); % print
%%%%%    QZ  = 1;  %  ie IGNORE the above QZ


        % The syy1,syy2,phi have compressed values u1,u2,u3 in the
        % [0,1]X[0,1]X[0,1] cube. One then uses
        % trilinear interpolation on the TBLKPHI table to get
        % K(syy1,syy2,phi)

        %  lookupKphi(phi); % gives the BACKPROJECTION of this kernel.
        %  Recall that ff starts out as zero
%-----------------------------
          % Thus center of scanner is at Cz = -LLL+DZ*count; and scanner detects all z in [Cz-LL,Cz+LL]
          % and we want z1,z2 inside. That is, max([abs(z1-Cz),abs(z2-Cz)])<LL
%===============================================

%%         Cz = 0.0; % We do NOT translate thescanner  %Jan2_19
   Cz = -LLL+DZ*count;
     if (  max(  [abs(z1-Cz),abs(z2-Cz)]  )<LL  )
       f1t = f1t + QZ*lookup1taug14(phi); % final image is the sum of these backprojections
%      f1t = f1t + QZ*lookup1tfeb6(phi); % final image is the sum of these backprojections

       else     % do nothing; z1,z2 not inside scanner
       rejecLL = rejecLL+1;
       endif  % end of the LL if
             endif   % end of the D<0 if

  count = count + 1;
           end    % end of ksrc loop over the Nsrc emission sources
    Cz                 % print
    rejecLL = rejecLL  % print
    rejecD = rejecD    % print
    accepLL = Nsrc-rejecLL    % print
%   accepLL = Nmax-rejecLL    % print
    count = count     % print
    cpurecon = etime(clock,ttt1) % print

%=====================================================

      save svff.mat f1t PHIBMAX Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc cpurecon  aa bb cc Nsrc RR 


          clear all
          close all
          load svff.mat

  ff1t(:,:) = f1t(:,Nyy+1,:);
  mmax = max(max(ff1t));
  invmax = 1/mmax;
  ff1t = invmax*ff1t;


%%%%====================================================
     % Now some graphical output to visualize the reconstructed ff1 etc


%-----------------------
    hold on

  figure (2); plot(ff1t)
  colormap(jet); legend('ff1t');
    hold on

  figure (3); plot(min(.2,ff1t));
  colormap(jet); legend('ff1t');
    hold on

  figure (4); mesh(min(.2,ff1t));
  colormap(jet); legend('ff1t');
    hold on


%    END

