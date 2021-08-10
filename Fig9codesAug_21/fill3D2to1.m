
            %  file fill3D2to1.m  % Sept1_19  

      function CC = fill3D2to1(BB);

  % We have a 3D matrix BB(j,k,p) with CORRECT image values on the odd j,k,p vertices.
  %     BB is Nx+1byNy+1byNz+1 with NX and Ny and Nz even.
  %     We want to fill in all values by bilinear interpolation.
  %     This takes from a dxx=2 to a dxx=1 image CC.

       sze = size(BB);
       Nx = sze(1);Ny = sze(2);Nz = sze(3);

        CC= BB;

%-----------------------------------------
  % 1st: fill in the odd xy-slices 
                  for p = 1:2:Nz+1
      CC(:,:,p) = fill2D2to1(BB(:,:,p));
                  end    

  % 2nd: fill in the even xy-slices by interpolation from the odd ones
                  for p = 2:2:Nz
       CC(:,:,p) = .5*CC(:,:,p-1) + .5*CC(:,:,p+1);
                  end
  % DONE!
%-----------------------------------------



  



