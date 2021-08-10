
   % file view124July24_21.m  (was view124Jan7_21.m

%%   save svff.mat f1 f1t f124 f12t f24tt f4 PHIBMAX Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc cpurecon  aa bb cc Nsrc  R12 R24 RG12 RG24 RTOF fac ;  % output string saved in svff.mat
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

  figure (2); plot(min(0.1,ff1));
  colormap(jet); legend('f1');
    hold on

  figure (3); plot(min(0.1,ff1t));
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

