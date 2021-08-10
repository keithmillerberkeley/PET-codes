
       % file view2Df124jan29_21.m  (was ~~ view2Df124mar31_20.m

          clear all
          close all
          load svff.mat

%  %%%%===================================================

%%    save svff.mat f1 f124 f12t f24tt f4 f2 PHIBMAX Ny1 Ny2 Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc cpurecon  aa bb cc Nsrc R1t R12 R24 RG12 RG24 RTOF N1 N2 N3 N4 N1234 shift1 shift2 shift3 shift4;       % output string saved in svff.mat

  PHIBMAX = 'pi/6'
  RTOF  % print
  R1t = 10
  R12  % print
  R24  % print
  fac
  RG12
  RG24
  Nsrc
  cpurecon
 
  ff1(:,:) = f1(:,:,Nzz+1);


  ff124(:,:) = f124(:,:,Nzz+1);

  ff12t(:,:) = f12t(:,:,Nzz+1);

  ff24tt(:,:) = f24tt(:,:,Nzz+1);

   ff2(:,:) = f2(:,:,Nzz+1);

   ff4(:,:) = f4(:,:,Nzz+1);

%%%%%%%%   we fill in from dxx = 2 or 4 to dxx =1 by 2D interpolation

  ff24ttfill = fill2D2to1(ff24tt);

  ff4fill = fill2D4to1(ff4);

  ff124fill = ff12t + ff24ttfill + ff4fill;

%%%%====================================================
     % Now some graphical output to visualize the reconstructed ff1 etc

  max1 = 1.030*max(max(ff1))
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

%-----------------------
    hold on

  figure (2); mesh(ff1);
  colormap(jet); legend('f1');
    hold on

  figure (3); mesh(ff2);
  colormap(jet); legend('f2');
    hold on

  figure (4); mesh(ff124fill);
  colormap(jet); legend('f124fill');
    hold on

  figure (5); mesh(min(.50,ff124fill));
  colormap(jet); legend('f124fill');
    hold on

  figure (6); plot(min(.20,ff124fill));
  colormap(jet); legend('f124fill');
    hold on

  figure (7); mesh(ff12t);
  colormap(jet); legend('f12t');
    hold on

     hold off

