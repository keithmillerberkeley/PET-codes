
   % file cnteventJuly18_21.m ( from ~/rrrandPETfig8_Jun30/

   close all
   clear all

%    save svff.mat f1 PHIBMAX Ny1 Ny2 Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc Nmax aa bb cc TOFu N1 N2 N3 N4 N1234 cpurecon mmmaxf1  % output string saved in svff.mat  by randPET1

   load svff.mat                % Nxx = 38 38 05  Nsrc=1.0e3
   Nsrc = Nsrc    % print
   avINTPRF = 28.350    %  from July3_21 Notes 
   szef1 = size(f1)   % print  77X77X11

   integral = dxx*dyy*dzz*sum(sum(sum(f1)))
   avintegral = integral/(1.6*Nsrc)  % print The 1.6 is needed july18_21
   integralneg = dxx*dyy*dzz*sum(sum(sum(min(0,f1))))
   intnegOint = integralneg/integral

   mmmaxf1 = max(max(max(f1)))  % print
   mmminf1 = min(min(min(f1)))   % print July18_21 better
   minOmaxf1 = mmminf1/mmmaxf1   % print
   ff1(:,:) = f1(:,:,Nzz+1);
   ff1(:,:) = (1/mmmaxf1)*ff1;


  colormap(jet)
  figure(2); contour(ff1,20);colormap(jet); colorbar
  avINTf1 = sum(sum(sum(f1)))/(1.6*Nsrc)   % print
  aavINTf1 = avINTf1/avINTPRF  % print

  BB = f1(49:69,15:35,:);
   szeBB = size(BB)
  avINTBB = sum(sum(sum(BB)))/(0.5*Nsrc)   % print
  aavINTBB = avINTBB/avINTPRF  % print
  BBB = BB(:,:,Nzz+1);
  figure(3); contour(BBB,20);colormap(jet); colorbar
  
  CC = f1(41:61,41:61,:);
   szeCC = size(CC)
  avINTCC = sum(sum(sum(CC)))/(1.0*Nsrc)   % print
  aavINTCC = avINTCC/avINTPRF   % print
  CCC = CC(:,:,Nzz+1);
  figure(4); contour(CCC,20);colormap(jet); colorbar
  
  DD = f1(15:35,49:69,:);
   szeDD = size(DD)
  avINTDD = sum(sum(sum(DD)))/(0.1*Nsrc)   % print
  aavINTDD = avINTDD/avINTPRF  % print
  DDD = DD(:,:,Nzz+1);
  figure(5); contour(DDD,20);colormap(jet); colorbar
  
  EE = f1(1:37,1:37,:);
   szeEE = size(EE)
  avINTEE = sum(sum(sum(EE)))/(1.0*Nsrc)   % print
  aavINTEE = avINTEE/avINTPRF   % print
  EEE = EE(:,:,Nzz+1);
  figure(6); contour(EEE,20);colormap(jet); colorbar
  
  figure(7); mesh(ff1);colormap(jet); 


 
