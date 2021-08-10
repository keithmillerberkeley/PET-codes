
            %  file fill3D4to1.m  % Sept2_19  

      function CC = fill3D4to1(BB);

  % We have a 3D matrix BB(j,k,p) with CORRECT image values on 1,5,9,... vertices.
  %     BB is Nx+1byNy+1byNz+1 with NX and Ny and Nz multiples of 4.
  %     We want to fill in all values by bilinear interpolation.
  %     This takes from a dxx=4 to a dxx=1 image CC.

       sze = size(BB);
       Nx = sze(1);Ny = sze(2);Nz = sze(3);

        CC= BB;

%-----------------------------------------
  % 1st: fill in the p=1,5,9,... xy-slices 
                  for p = 1:4:Nz+1
   CC(:,:,p) = fill2D4to1(BB(:,:,p));
                  end
%-----------------------------------------------------
% 2nd: fill in the p=2,3,4;6,7,8; ... xy-slices by interp from the 1,5,9,... ones
                for p = 1:4:Nz-3
       CC(:,:,p+1) = .75*CC(:,:,p) + .25*CC(:,:,p+4);
       CC(:,:,p+2) = 0.5*CC(:,:,p) + 0.5*CC(:,:,p+4);
       CC(:,:,p+3) = .25*CC(:,:,p) + .75*CC(:,:,p+4);
                  end
  % DONE!
%-----------------------------------------



  



