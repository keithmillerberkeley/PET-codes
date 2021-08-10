
  % file view124July27_21.m  (was view124Jan20_21.m

     clear all
     close all
   load svff.mat

%%   save svff.mat f124 f12t f24tt f4 f2 PHIBMAX Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc cpurecon  aa bb cc Nsrc  R12 R24 RG12 RG24 TOFu fac ;  % output string saved in svff.mat

  PHIBMAX
  Nsrc
  R12 
  R24  
  TOFu  
  fac
  aa  
  bb  
  cc
  cpurecon

   ff12t(:,:) = f12t(:,Nyy+1,:);
  ff24tt(:,:) = f24tt(:,Nyy+1,:);
   ff4(:,:) = f4(:,Nyy+1,:);
  ff124(:,:) = f124(:,Nyy+1,:);

   ff2(:,:) = f2(:,Nyy+1,:);
   mmax124 = 2.483e4;
   mmax2 = 7268;
   ratio = mmax2/mmax124;
   ff2 = ratio*ff2;
  
  ff24ttfill = fill2D2to1(ff24tt);
  ff4fill = fill2D4to1(ff4);

  ff124fill = ff12t + fac*ff24ttfill + fac*fac*ff4fill;


%%%%====================================================
     % Now some graphical output to visualize the reconstructed ff1 etc


%-----------------------
    hold on

  figure (2); plot(ff12t);
  colormap(jet); legend('ff12t');
    hold on

  figure (3); plot(ff24tt);
  colormap(jet); legend('ff24tt');
    hold on

 figure (4); plot(ff4);
  colormap(jet); legend('ff4');
    hold on

  figure (6); plot(ff124fill);
  colormap(jet); legend('ff124fill');
    hold on

  figure (7); plot(min(0.5,ff124fill));
  colormap(cubehelix); legend('ff124fill');
% colormap(jet); legend('ff124fill');
    hold on

  figure (8); plot(min(0.5,ff124));
  colormap(jet); legend('ff124');
    hold on

  figure (9); contour(ff124fill);
  colormap(jet); legend('ff124fill');
    hold on

  figure (10); mesh(ff124fill);
  colormap(jet); legend('ff124fill');
    hold on

  figure (11); mesh(ff124);
  colormap(jet); legend('ff124');
    hold on

  figure (13); mesh(ff2);
  colormap(jet); legend('ff2');
    hold on

  figure (12); mesh(min(.5,ff124fill));
  colormap(jet); legend('min(.5,ff124fill)');
    hold on

%    END
