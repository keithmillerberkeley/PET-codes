
     % in ~/newFig12feb11_21/NNnewFig12feb24_21/
   
         % file mk2X2newFig12july27_21.m (was mk2X2newFig12feb25_21.m  (was mk2X2newFig11feb6.m (was mk2X2newFig11.m  (was file mk2X2Fig10.m
 close all
 clear all

  
%--------------------------------

  load svff.matJuly27_21AD; % sin(phib) atan(15/40)R12=8 RTOF=10 fac=.8300 Nsrc=1.0e7 5cubes
  ff124(:,:) = f124(:,Nyy+1,:);   % This was all on dxx=1 grid; 53x1x53
  ff12t(:,:) = f12t(:,Nyy+1,:);   % This was all on dxx=1 grid
  ff24tt(:,:) = f24tt(:,Nyy+1,:); % This was all on dxx=1 grid
  ff4(:,:) = f4(:,Nyy+1,:);       % This was all on dxx=1 grid
  ffAD = ff12t + fac*fill2D2to1(ff24tt) + fac^2*fill2D4to1(ff4); %Good,
  mmax = max(max(ffAD)); immax = 1/mmax;
  ffAD = immax*ffAD;


  load svff.matJuly27_21B; % sin(phib) atan(15/40)R12=8 RTOF=10 fac=.8300 Nsrc=1.0e6 5cubes
  ff124(:,:) = f124(:,Nyy+1,:);   % This was all on dxx=1 grid; 53x1x53
  ff12t(:,:) = f12t(:,Nyy+1,:);   % This was all on dxx=1 grid
  ff24tt(:,:) = f24tt(:,Nyy+1,:); % This was all on dxx=1 grid
  ff4(:,:) = f4(:,Nyy+1,:);       % This was all on dxx=1 grid
  ffB = ff12t + fac*fill2D2to1(ff24tt) + fac^2*fill2D4to1(ff4); %Good,
  mmax = max(max(ffB)); immax = 1/mmax;
  ffB = immax*ffB;

  
  load svff.matJuly27_21C; % sin(phib) atan(15/40)R12=8 RTOF=10 fac=.8300 Nsrc=1.0e5 5cubes
  ff124(:,:) = f124(:,Nyy+1,:);   % This was all on dxx=1 grid; 53x1x53
  ff12t(:,:) = f12t(:,Nyy+1,:);   % This was all on dxx=1 grid
  ff24tt(:,:) = f24tt(:,Nyy+1,:); % This was all on dxx=1 grid
  ff4(:,:) = f4(:,Nyy+1,:);       % This was all on dxx=1 grid
  ffC = ff12t + fac*fill2D2to1(ff24tt) + fac^2*fill2D4to1(ff4); %Good,
  mmax = max(max(ffC)); immax = 1/mmax;
  ffC = immax*ffC;


%%%%%%%%%%%%%%%%%


   hold on
  figure(2); colormap(jet);

  subplot(2,2,1); mesh(ffAD);
  subplot(2,2,2); mesh(ffB); 
  subplot(2,2,3); mesh(ffC); 
  subplot(2,2,4); mesh(min(0.5,ffAD));


   hold off



