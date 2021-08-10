

   % file view1tjuly6_21.m  (was view1july5_21.m

%%    save svff.mat f1t PHIBMAX R1t LG1t Ny1 Ny2 Nphi dxx dyy dzz Nxx Nyy Nzz Nsrc Nmax aa bb cc TOFu N1 N2 N3 N4 N1234 mmmaxf1t mmminf1t minOmaxf1t cpurecon avintegral  % output string saved in svff.mat
          clear all
          close all
          load svff.mat

  R1t
  TOFu
  LG1t
  Nsrc
  cpurecon

  ff1t(:,:) = f1t(:,:,Nzz+1);
  ff1t = (1/mmmaxf1t)*ff1t;

   
   mmmaxf1t = max(max(max(f1t)))    % print
   mmminf1t = min(min(min(f1t)))
   minOmaxf1t = mmminf1t/mmmaxf1t

   intpos = sum(sum(sum(max(0.0,f1t))))
   intneg = sum(sum(sum(min(0.0,f1t))))
   int    = sum(sum(sum(f1t)))
   intnegOint = intneg/int


%%%%====================================================
     % Now some graphical output to visualize the reconstructed ff1 etc


%-----------------------
    hold on

  figure (2); mesh(ff1t);
  colormap(jet); legend('ff1t');
    hold on

  figure (3); plot(ff1t);
  colormap(jet); legend('ff1t');
    hold on

  figure (4); plot(min(0.2,ff1t));
  colormap(jet); legend('ff1t');
    hold on

  figure (5); contour(ff1t);
  colormap(jet); colorbar; legend('ff1t');
    hold on

  figure (6); image(60*ff1t);
  colormap(jet); colorbar;
% colormap(jet); colorbar; legend('ff1t');
    hold on

    hold off


%    END of the reconstruction of the image ff on 4 cubes of recon pts




