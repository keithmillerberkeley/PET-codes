
   % file cnteventJuly17_21.m  ( from ~/rrrandPETfig8_Jun30/

   close all
   clear all

   load svff.mat  % 
   R1t
   TOFu
   LG1t
   cpurecon
   Nsrc = Nsrc    % print
   avINTPRF = 15.000   % print  % calibrated Aug3_21 with dxx=.5 Nxx=24 24 24 Nsrc=1.0e4
%  avINTPRF = 13.217   % print  % calibrated July17_21

    szef1t = size(f1t)   % print  79X79X11
  
  ff1t = f1t(:,:,Nzz+1);
  colormap(jet)
  figure(2); contour(ff1t);colormap(jet); colorbar
  avINTf1t = sum(sum(sum(f1t)))/(1.6*Nsrc)   % print
  aavINTf1t = avINTf1t/avINTPRF  % print

  BB = f1t(49:69,15:35,:);
   szeBB = size(BB)
  avINTBB = sum(sum(sum(BB)))/(0.5*Nsrc)   % print
  aavINTBB = avINTBB/avINTPRF  % print
  BBB = BB(:,:,Nzz+1);
  figure(3); contour(BBB,20);colormap(jet); colorbar
  
  CC = f1t(41:61,41:61,:);
   szeCC = size(CC)
  avINTCC = sum(sum(sum(CC)))/(1.0*Nsrc)   % print
  aavINTCC = avINTCC/avINTPRF   % print
  CCC = CC(:,:,Nzz+1);
  figure(4); contour(CCC,20);colormap(jet); colorbar
  
  DD = f1t(15:35,49:69,:);
   szeDD = size(DD)
  avINTDD = sum(sum(sum(DD)))/(0.1*Nsrc)   % print
  aavINTDD = avINTDD/avINTPRF   % print
  DDD = DD(:,:,Nzz+1);
  figure(5); contour(DDD,20);colormap(jet); colorbar
  
  EE = f1t(1:39,1:39,:);
   szeEE = size(EE)
  avINTEE = sum(sum(sum(EE)))/(1.0*Nsrc)   % print
  aavINTEE = avINTEE/avINTPRF   % print
  EEE = EE(:,:,Nzz+1);
  figure(6); contour(EEE,20);colormap(jet); colorbar
  


 
