
       %  file mk2X2newFig6July6_21.m  (was mk2X2newFig10Feb2_21.m  (was mk2X2F1g13Aug28.m

  close all
  clear all
%------------------------------

% load svff.matJuly5_21D % this was 49X49X49 so I deleted it Aug4_21 
% ff1(:,:) = f1(:,Nyy+1,:);
% ff1p = ff1'; %this puts z across in "plot"
  load svff1pD.mat    % this is 49X49 version

  subplot(2,2,1); 
  plot(ff1pD);
% plot(ff1p);
% contour(ff1p,20);
  colormap(jet);
    hold on
 
%--------------------------

% load svff.matJuly5_21A % this was 49X49X49 so I deleted it Aug4_21
% load svff.matJuly5_21A  % without the width-factor
% ff1(:,:) = f1(:,Nyy+1,:);
% ff1p = ff1'; %this puts z across in "plot"
  load svff1pA.mat    % this is 49X49 version

  subplot(2,2,2); 
  plot(ff1pA);
% plot(ff1p);
% contour(ff1p,20);
  colormap(jet);
    hold on

%--------------------------

  load svff.matJuly21_21B
  fff1(:,:) = f1(:,:,Nzz+1);
  fff1 = (1/mmmaxf1)*fff1;
  

  subplot(2,2,3); 
  mesh(fff1);
% mesh(min(.5,fff1));
  colormap(jet);
    hold on

%--------------------------

  subplot(2,2,4); 
  plot(min(.20,fff1));
  colormap(jet);
    hold on
  
%--------------------------
    hold off
