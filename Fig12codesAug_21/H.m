
            %       file H.m

   function  Ht = H(t)

%--------------------------------------------------------
% evaluate the 2-D convolution kernel H(t) by its Fourier integral for small t,
% say t < 7


%  if (abs(t) < 7)
   if (abs(t) <30)
%  if (abs(t) <10)

%rmax = 7;
 rmax = 30;
%nrper = 9*7;
 nrper = 29*7;

 dr = 1/nrper;
 nr = nrper*rmax+1;
 r = [0:dr:rmax];

 funcr = zeros(1,nr);
funcr = r.*exp(-.5*r.*r).*cos(t*r);  % nr vector
%    for uu = 1:nr
%rr = r(uu);
%funcr(uu) = rr *exp(-.5*rr*rr)*cos(t*rr);
%funcr(uu) = rr*HC(rr)*cos(t*rr);
%    end   %  end of uu loop

 funcr(1) = .5*funcr(1);    % .5 at endpt; it's already zero
 funcr(nr) = .5*funcr(nr);  % .5 At endpt; it's already almost zero

 Ht = dr*sum(funcr) ;         % integral by Trapezoid Rule

% Aug 21, 2004 I add the Euler-MacLaurin end-pt correction for the
%Trap rule. It is, for this integral on [0,inf) of f = funcr
% correction = f'*dr^2/(6*2!) + f'''*dr^4/(-30*4!) + f'''''*dr^6/(42*6!) + ...

%-------------------------------------------------------------------
    else
% asymptotic expansion, from successive integ by parts, for large t

   Ht = -1*t^-2  - 3*t^-4  -15*t^-6  -105*t^-8 ;
%  Ht = -1*t^-2  - 3*t^-4  -15*t^-6  -105*t^-8 ;
%%%%    Ht = -1*t^-2  ;

%-------------------------------------------------------
    end



