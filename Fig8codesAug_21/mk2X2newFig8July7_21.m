
       %  file mk2X2newFig8July7_21.m  (was mk2X2newFig10Feb2_21.m  (was mk2X2F1g13Aug28.m

  close all
  clear all
%------------------------------

  load svff.matAug2_21A  % Nsrc=1.0e6 pi/6 R1t=20 TOFu=10 Truncates at R1t
% load svff.matJuly19_21A  % Nsrc=1.0e6 pi/6 R1t=20 ERROR, does NOT truncate at R1t
  
  R1t
  TOFu
  LG1t
  cpurecon
  Nsrc


  ff1t(:,:) = f1t(:,:,Nzz+1);
  mmmaxf1t = max(max(max(f1t)))
  ff1t = (1/mmmaxf1t)*ff1t;

  subplot(2,2,1); 
  mesh(min(.5,ff1t));
  colormap(jet);
    hold on
 
%--------------------------

  subplot(2,2,2); 
  contour(ff1t,20); colorbar;
  colormap(jet);
    hold on

%--------------------------

  subplot(2,2,3); 
 plot(min(.2,ff1t));
  colormap(jet);
    hold on

%--------------------------

  load svff.matAug2_21B    % Nsrc=1.0e6 R1t=10 TOFu=20 positive bulge
% load svff.matJuly19_21B  % Nsrc=1.0e4 R1t=20 TOFu=10 

  ff1t(:,:) = f1t(:,:,Nzz+1);
  mmmaxf1t = max(max(max(f1t)));
  ff1t = (1/mmmaxf1t)*ff1t;

  subplot(2,2,4); 
  plot(min(.2,ff1t));
% mesh(min(.50,ff1t));
  colormap(jet);
    hold on
  
%--------------------------
    hold off
