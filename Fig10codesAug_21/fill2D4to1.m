
            %  file fill2D4to1.m  % Sept1_19  

      function CC = fill2D4to1(BB);

  % We have a 2D matrix BB(j,k) with CORRECT image values on the 1,5,9,... j,k vertices.
  %     BB is Nx+1byNy+1 with Nx and Ny multiples of 4.
  %     We want to fill in all values by bilinear interpolation.
  %     This takes from a dxx=4 to a dxx=1 image CC.

       sze = size(BB);
       Nx = sze(1);Ny = sze(2);

        CC= BB;
  % 1st: fill in the 1,5,9,... x-lines.
                   for k = 1:4:Ny+1
                for j = 1:4:Nx-3
       CC(j+1,k) = .75*BB(j,k) + .25*BB(j+4,k);
       CC(j+2,k) = .5*BB(j,k) + .5*BB(j+4,k);
       CC(j+3,k) = .25*BB(j,k) + .75*BB(j+4,k);
                end;
                   end;
  % 2nd: fill in the 2,3,4;6,7,8; ... x-lines
                for k = 1:4:Ny-3
       CC(:,k+1) = .75*CC(:,k) + .25*CC(:,k+4);
       CC(:,k+2) = .5*CC(:,k) + .5*CC(:,k+4);
       CC(:,k+3) = .25*CC(:,k) + .75*CC(:,k+4);
                  end
  % DONE!




  



