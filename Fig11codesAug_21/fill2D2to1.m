
            %  file fill2D2to1.m  % Sept1_19  

      function CC = fill2D2to1(BB);

  % We have a 2D matrix BB(j,k) with CORRECT image values on the odd j,k vertices.
  %     BB is Nx+1byNy+1 with Nx and Ny even.
  %     We want to fill in all values by bilinear interpolation.
  %     This takes from a dxx=2 to a dxx=1 image CC.

       sze = size(BB);
       Nx = sze(1);Ny = sze(2);

        CC= BB;
  % 1st: fill in the odd x-lines.
                   for k = 1:2:Ny+1
                for j = 2:2:Nx
       CC(j,k) = .5*BB(j-1,k) + .5*BB(j+1,k);
                end;
                   end;
  % 2nd: fill in the even x-lines. 
                  for k = 2:2:Ny
       CC(:,k) = .5*CC(:,k-1) + .5*CC(:,k+1);
                  end
  % DONE!




  



