
    %  mk2X2newFig11july26.m  (was mk2X2newFig11feb6.m (was mk2X2newFig11.m  (was file mk2X2Fig10.m
 close all
 clear all

  
%--------------------------------

  load svff.matJuly26AB % sin(phib) atan(.8*15/40)R12=8 RTOF=10 fac=.8200 Nsrc=1.0e7 aa=2 2 45
% load svff.matJan22_21B; % sin(phib) atan(.8*15/40)R12=8 RTOF=10 fac=.82097 Nsrc=1.0e7 aa=2 2 45
  ff124(:,:) = f124(:,Nyy+1,:);   % This was all on dxx=1 grid
  ff12t(:,:) = f12t(:,Nyy+1,:);   % This was all on dxx=1 grid
  ff24tt(:,:) = f24tt(:,Nyy+1,:); % This was all on dxx=1 grid
  ff4(:,:) = f4(:,Nyy+1,:);       % This was all on dxx=1 grid
  ffAB = ff12t + fac*fill2D2to1(ff24tt) + fac^2*fill2D4to1(ff4); 
  mmax = max(max(ffAB)); immax = 1/mmax;
  ffAB = immax*ffAB;


  
  load svff.matJuly26C; % sin(phib) atan(15/40) R12=8 RTOF=10 fac=.8300 Nsrc=1.0e7 aa=2 2 45
% load svff.matJan21_21A; % sin(phib) atan(15/40) R12=8 RTOF=10 fac=.8269 Nsrc=1.0e7 aa=2 2 45
  ff12t(:,:) = f12t(:,Nyy+1,:);
  ff24tt(:,:) = f24tt(:,Nyy+1,:);
  ff4(:,:) = f4(:,Nyy+1,:);
  ffC = ff12t + fac*fill2D2to1(ff24tt) + fac^2*fill2D4to1(ff4);
  mmax = max(max(ffC)); immax = 1/mmax;
  ffC = immax*ffC;
  

  load svff.matJuly26D; % sin(phib) atan(.8*15/40) R1t=20 RTOF=10 Nsrc=1.0e7 aa=2 2 45
% load svff.matFeb6_21A; % sin(phib) atan(.8*15/40)R1t=20 RTOF=10 Nsrc=1.0e7 aa=2 2 45
  ffD(:,:) = f1t(:,Nyy+1,:);   % This was all on dxx=1 grid
  mmax = max(max(ffD)); immax = 1/mmax;
  ffD = immax*ffD;
  mmax = max(max(ffD)); immax = 1/mmax;
  ffD = immax*ffD;



   hold on
  figure(2); colormap(jet);

  subplot(2,2,1); contour(ffAB);legend('f124fill');
  subplot(2,2,2); plot(ffAB);
  subplot(2,2,3); plot(ffC);
  subplot(2,2,4); plot(ffD);

   hold on
  figure(5); colormap(jet);
  mesh(min(.1,ffAB)); legend('f124fill atan(.8*15/40)');

   hold on
  figure(6); colormap(jet);
  mesh(min(.1,ffC));legend('f124 atan(15/40)');

   hold on
  figure(7); colormap(jet);
  mesh(min(.1,ffD));legend('f1t atan(.8*15/40)');


   hold off



