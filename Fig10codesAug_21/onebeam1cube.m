
          %  file onebeam1cube.m 
        % this produces a single PET source beam with abs(phi) < PHIBMAX.

  global PHIBMAX;  % max phi (angle from equator) of our kernel
  global aa  bb  cc;
  global ssrc  dir;
   
  global count     % sept2_17  we keep track of number of src's
  global TOFu      % TOF Gaussian uncertainty

  sx = -aa + 2*aa*rand; sy = -bb + 2*bb*rand; sz = -cc+ 2*cc*rand;
  src = [sx,sy,sz];  % the position vector

  dirx = -1+2*rand;
  diry = -1+2*rand;
  dirz = -1+2*rand;
  dir =  [dirx,diry,dirz]; % random points in the [-1,+1] cube
  rr = norm(dir);
  dir = dir/rr;
  dir = sign(.00000123+dir(3))* dir;   % unit direction 3-vector with

     while ( (rr < .1) || (rr > 1) || (abs(dir(3)) > .9999999*sin(PHIBMAX)) )
              %   NOTE the "||" is an "or"
  dirx = -1+2*rand;
  diry = -1+2*rand;
  dirz = -1+2*rand;
  dir =  [dirx,diry,dirz];
  rr = norm(dir);
  dir = dir/rr; % the unit direction vector
      end  % end of the while

%  Aug11_17 I shift the source pt to another pt on the beam (so that we don't
%    know WHERE on the beam the source is).   ssrc is the "shifted source".
  ssrc = src + TOFu*randn*dir;  %  goes with  random normal 

