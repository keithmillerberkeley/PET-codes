
          %  file onebeam5cubeJuly27_21.m  (was onebeam5cubeOct4.m  from Mar20.mfrom Jan6.m
        % this produces a single PET source beam with abs(phi) < PHIBMAX.

  global PHIBMAX;  % max phi (angle from equator) of our kernel
  global aa  bb  cc;
  global ssrc  dir;
  global ssrc  dir;
  global count     % sept2_17  we keep track of number of src's

  global TOFu R12 R24 RG12 RG24
  global shift1 shift2 shift3 shift4 shift5;
  global N1 N2 N3 N4 N5 N12345;


% we get a dir with abs(phi) < PHIBMAX 
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
%%%%%%%%%%%%%%

%  We randomly put sources in one of five cubes NOTE aa = 2 2 2 Mar3_21
  sx = -aa + 2*aa*rand; sy = -bb + 2*bb*rand; sz = -cc + 2*cc*rand;

  sssrc = [sx,sy,sz];  % the position vector

  Rand = N12345*rand;
      if (Rand <= N1);
   src = sssrc + shift1;
 
      elseif and(Rand >N1 , Rand <= (N1+N2));
   src = sssrc + shift2;
     
      elseif and(Rand >(N1+N2) , Rand <= (N1+N2+N3));
   src = sssrc + shift3;

      elseif and(Rand >(N1+N2+N3) , Rand <= (N1+N2+N3+N4));
   src = sssrc + shift4;

      elseif and(Rand >(N1+N2+N3+N4) , Rand <= (N1+N2+N3+N4+N5));
     sssrc = [7.5*sx 2*sy 17.5*sz]; % temporary test Feb19_21
%    sssrc = [20*sx 2*sy 20*sz];
   src = sssrc + shift5;

%     else

      endif
     
%%%%%%%%%%%%%%

  ssrc = src + TOFu*randn*dir;

