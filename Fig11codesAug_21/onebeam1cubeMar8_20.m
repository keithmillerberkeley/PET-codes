
          %  file onebeam1cubeMar8_20.m  , from onebeam1cubeDec9_19.m
        % this produces a single PET source beam with abs(phi) < PHIBMAX.

  global PHIBMAX;  % max phi (angle from equator) of our kernel
  global aa  bb  cc;  % from randPET3D
  global ssrc  dir; % generated here; passes to vecmap
%  global shift; % from randPET3D
   
  global count RTOF     % sept2_17  we keep track of number of src's

%        sx = -aa+ 2*aa*rand;% sy = -bb + 2*bb*rand;% sz = -cc + 2*cc*rand;

       ssx = -aa+ 2*aa*rand; sy = -bb + 2*bb*rand; ssz = -cc + 2*cc*rand;
       sx = .70711*ssx + .70711*ssz; sz = -.70711*ssx + .70711*ssz;

%      ssx = -aa+ 2*aa*rand; sy = -bb + 2*bb*rand; ssz = -cc + 2*cc*rand;
%      sx =  .86603*ssx + .50000*ssz; sz = -.50000*ssx + .86603*ssz;

  src = [sx,sy,sz] ;  % the position vector

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
%ssrc = src + randn*dir;   % Dec9,2019
 ssrc = src + TOFu*randn*dir; %Dec29,2019 this is the "presumed source psrc"

