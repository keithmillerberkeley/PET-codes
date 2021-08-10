
         % mk2X2newFig9jan31_21.m  (was mk2X2Fig12.m  Apr15_20
 close all
 clear all

%--------------------------------
 load svTBLK124oct17  % works. Found in ~/savePETTOFjul16_19/mgrid124Oct11_19
%load svTBLK124apr13_20  % was ~~svTBLK124oct17

  figure(2); colormap(jet);
  subplot(2,2,1);
     plot(00+rrr3.*TBLK1(:,:,01) );
         hold on
     plot(00+rrr3.*TBLK2(:,:,01) );
         hold on
     plot(00+rrr3.*TBLK4(:,:,01) );
         hold on

     plot(05+rrr3.*TBLK1(:,:,03) );
         hold on
     plot(05+rrr3.*TBLK2(:,:,03) );
         hold on
     plot(05+rrr3.*TBLK4(:,:,03) );
         hold on

     plot(10+rrr3.*TBLK1(:,:,04) );
         hold on
     plot(10+rrr3.*TBLK2(:,:,04) );
         hold on
     plot(10+rrr3.*TBLK4(:,:,04) );
         hold on
   legend('(K1,K2,K4)*r**3  dr=1', 'phi = (.00 .50 .75)*pi/6');
%  legend('(K1,K2,K4)*rrr3  dr=1', 'phi = (.00 .50 .75)*pi/6');
         hold on
%-------------------------------------
 load svff.matJuly23_21B %4cubes Nsrc=1.0e6 R1t=20 R12=8 TOFu=10 tuned fac=.8200
%load svff.mat  % 4cubes  Nsrc=1.0e3
%%%%% load svff.matApr5A  % 4cubes  Nsrc=1.0e7 R1t=10 R2t=20 TOFu=10
%%%%% load svff.matApr6C  % 4cubes rotated Nsrc=1.0e7
% load svff.matJan15A_21  % 4cubes  Nsrc=1.0e6 R12=5 TOFu=10 fac=.7141
%% load svff.matApr5B  % 4cubes  Nsrc=1.0e6 R1t=10 R2t=20 TOFu=10
%%%%% load svff.matApr1B      % 4cubes  Nsrc=1.0e6 R1t=20 R2t=40 TOFu=10
  max1 = 1.01*max(max(f1));
  inv1 = 1/max1;

  ff24ttfill = fill2D2to1(f24tt);
  ff4fill = fill2D4to1(f4);

  ff24ttfill = inv1*ff24ttfill;
  ff4fill = inv1*ff4fill;
  ff12t = inv1*f12t;

  ff124fill = ff12t + ff24ttfill + ff4fill;
%---------------------------------
 
 figure(2); colormap(jet);
  subplot(2,2,2);
%%%%%   mesh(ff124fill);
   mesh(min(0.5,ff124fill));
    hold on

%-----------------------------------
 figure(2); colormap(jet);
  subplot(2,2,4);
   mesh(ff12t);
    hold on

%-----------------------------------
 figure(2); colormap(jet);
  subplot(2,2,3);
   plot(min(0.2,ff124fill));
    hold on


   hold off



