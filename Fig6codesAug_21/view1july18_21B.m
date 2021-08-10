
%    file view1july18_21B.m   

%    save svff.mat f1 PHIBMAX Ny1 Ny2 Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc Nmax aa bb cc TOFu N1 N2 N3 N4 N1234 cpurecon   % output string saved in svff.mat  by randPET1

    clear all
    close all

   load svff.mat
   PHIBMAX
   TOFu
   Nsrc = Nsrc    % print
   cpurecon
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




%%%%====================================================
     % Now some graphical output to visualize the reconstructed ff1 etc

%-----------------------
    hold on

  figure (2); mesh(ff1);
  colormap(jet); legend('ff1');
    hold on

  figure (3); plot(ff1);
  colormap(jet); legend('ff1');
    hold on

  figure (4); plot(min(0.2,ff1));
  colormap(jet); legend('ff1');
    hold on

  figure (5); contour(ff1);
  colormap(jet); colorbar; legend('ff1');
    hold on

  figure (6); image(60*ff1);
  colormap(jet); colorbar;
    hold on

  figure (7); mesh(min(.5,ff1));
  colormap(jet); legend('ff1');
    hold on

%    END of the reconstruction of the image ff on 4 cubes of recon pts


