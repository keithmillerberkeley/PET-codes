

          %   file vecmapsep24.m 

%   This is called from randPET3D.m It computes the OUTPROJECTED screen
    %  coords syy1,syy2 from the xx,yy,zz reconstruction mesh, for
    %  the single random beam with emission at pt src in direction dir

    %  But in onebeam.m we shifted to the point ssrc = src+3*dir along
    %  the beam to hide the exact location of ssrc on the beam

 function  vecmapsep24(ssrc,dir)
 
  global  xx yy zz  ; % I pass these BIG input arrays as globals
  global  syy1 syy2 ; % I pass these BIG output arrays as globals
                    %  src is the   xyz position of the SINGLE source
                    %  dir is the unit direction vector of our SINGLE beam
  global  syy0;     % Jul11_18  added for TOF
  global  rrsq      % Sep24_19  added for TOF

  ww = dir;   % unit beam-direction 3-vector with positve dir(3)
  ww2 = [0,0,1];  % this points toward the NORTH POLE
  ww1 = [-1.153794071,1,0]; % Mar17,07   I need
      % to make this slightly random to ensure that ww,ww2,ww1 are not
      % accidentally linearly dependent

  %  We orthogonalize ww,ww2,ww1  by Gram-Schmidt
  ww2 = ww2 - (ww2*ww')*ww;
  ww2 = ww2/norm(ww2);

  ww1 = ww1 - (ww1*ww')*ww -(ww1*ww2')*ww2;
  ww1 = ww1/norm(ww1);


%------------------------------------------------
  % now ww2,ww1 are orthogonal to the beam, and tangent to the celestial sphere 
         % in the directions of increasing phi and theta. 

%%-----------------------------------------------

% Sep27-2016 I translate ssrc to be the origin
  utxx = xx - ssrc(1); % Feb24 x-components of vectors reconpts - ssrc
                            % NOTE matrix-scalar  has size of the matix
  utyy = yy - ssrc(2); % Feb24 y-components
  utzz = zz - ssrc(3); % Feb24 z-components

  syy1 = ww1(1)*utxx + ww1(2)*utyy + ww1(3)*utzz;
  syy2 = ww2(1)*utxx + ww2(2)*utyy + ww2(3)*utzz;

  syy0 = ww(1)*utxx + ww(2)*utyy + ww(3)*utzz;  % Jul11_18 addedfor TOF
  rrsq  = syy0.^2 +syy1.^2 +syy2.^2;   % Sep24_19  added for TOF

%  NOTE syy1 is coord in direction of increasing theta; across on the K screen 
%  NOTE syy2 is coord in direction of increasing phi; up on the K screen 

%  NOTE syy0 is coord in direction of beam 


%%% Matrices syy1 syy2 syy0 now returned as global variables; sep24 now also rrsq

