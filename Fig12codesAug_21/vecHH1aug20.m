
         %  file vecHH1aug20.m

%% This is a code for finding the  beta-rotated-HH(yy1,yy2,phi,phib) kernel
%%    on the yy1,yy2 screen  for given phi and phib
%%-----------------------------------------------------------------


       function HH = vecHH1aug20(phi,phib)

 global PHIBMAX;  % max phi (angle from equator) of our kernel
 global TBLHu;  % table of Hu(u) values, done once in main code 
 global NHu;      % # of TBLHu entries, usually 2001
 global yy1;      
 global yy2;   
 global Ny1;    
 global Ny2; 
% HH  this will contain the kernel on the yy1,yy2 screen for backprojection 

% Here yy1,yy2,HH are matrices of size Ny1 by Ny2

%---------------------------------
 sss = sqrt((sin(phib))^2 -(sin(phi))^2)/cos(phi);  %sine of rotation angle beta
 ccc = sqrt(1-sss*sss); % cosine of the rotation angle beta

% Note u,phib, sss and ccc are scalars

 invPMX = 1/PHIBMAX;  %  Nov3,07 avoid many divides

%=================================================
%Now compute the rotated kernel on the yy1,yy2 screen
%  tt,ss are the ROTATED coords of the matrix of y1,y2 points  

%   June23 I vectorize this with a table lookup for H(abs(t))
 
 tt=  ccc*yy1 + sss*yy2;
 ss = -sss*yy1 + ccc*yy2; % yy1,yy2,tt,ss are matrices 
 tt = abs(tt); % by symmetry in K we need only nonnegative tt,ss
 ss = abs(ss);
 
%tt = .5*tt;   % sep19-17 
%ss = .5*ss;   % MAY9_18
 
%tt = .25*tt;   % sep19-17 
%ss = .25*ss;
 
 % We do the linear interpolation of Hu = H(t) values from table lookup
 %    Recall t = 4u/(1-u) in making the 1xNHu TBLHu
 %    Thus u = t/(4+t) here; in [0,1)
  utt =  (NHu-1)*tt./(4+tt); % in [0,NHu-1)
  flutt = floor(utt);        % in [0,NHu-2]; thus 2+flutt is in [2,NHu]
  remutt = utt-flutt;


 Htt = (1-remutt).*TBLHu(1+flutt) + remutt.*TBLHu(2+flutt);
 Htt = Htt./(1+tt.*tt);% recall TBLHu gives H(t)*(1+t^2) instead of H(t)
 Hts = Htt.*(exp(-.5*ss.*ss)); %this is the beta_rotated HH on the screen

%----------------------------------------------------

%  now we rotate in the OPPOSITE angle and ADD the two
  tt=  ccc*yy1 - sss*yy2;  % yy1,yy2,tt,ss are matrices
 ss =  sss*yy1 + ccc*yy2;
 tt = abs(tt);
 ss = abs(ss);
 
%tt = .5*tt;   % sep19-17 
%ss = .5*ss;
 
%tt = .25*tt;   % sep19-17 
%ss = .25*ss;
 
 % We do the linear interpolation of Hu = H(t) values from table lookup
 %    Recall t = 4u/(1-u) in making the 1xNHu TBLHu
 %    Thus u = t/(4+t) here; in [0,1)
  utt =  (NHu-1)*tt./(4+tt);
  flutt = floor(utt);
  remutt = utt-flutt;


 Htt = (1-remutt).*TBLHu(1+flutt) + remutt.*TBLHu(2+flutt);
 Htt = Htt./(1+tt.*tt);% recall TBLHu gives H(t)*(1+t^2) instead of H(t)
 Hts2 = Htt.*(exp(-.5*ss.*ss)); %this is the -beta_rotated HH on the screen

    wphiphib = 1/sqrt((cos(phi))^2-(cos(phib))^2);
    HH = wphiphib*(Hts+Hts2);
%%     HH = (Hts+Hts2); %%% AUG29_17 I try WITHOUT the width factor
             % and with pi/3 I get a -.0206 undershoot of lip of peak

%========================================================
