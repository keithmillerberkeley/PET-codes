
       %  file mk2X2newFig10Feb2_21.m  (was mk2X2F1g13Aug28.m

  close all
  clear all
%------------------------------

  load svff.matJuly24_21A
% load svff.matJun29_21A
% load svff.matFeb1_21A
  f24ttfill = fill2D2to1(f24tt);
  f4fill = fill2D4to1(f4);
  f124fill = f12t + fac*f24ttfill + fac*fac*f4fill;

  subplot(2,2,1); 
   plot(min(.10,f1 )+.00);
  colormap(jet);% legend('f1');
    hold on
   plot(min(.10,f124fill )+0.15);
  colormap(jet);% legend('f124fill');
    hold on
   plot(min(.10,f1t )+0.30);
  colormap(jet); legend('f1 f124fill f1t');
    hold on
%--------------------------
  load svff.matJuly24_21B
% load svff.matJun29_21B
% load svff.matFeb1_21B
  f24ttfill = fill2D2to1(f24tt);
  f4fill = fill2D4to1(f4);
  f124fill = f12t + fac*f24ttfill + fac*fac*f4fill;

  subplot(2,2,2); 
   plot(min(.10,f1 )+.00);
  colormap(jet);% legend('f1');
    hold on
   plot(min(.10,f124fill )+0.15);
  colormap(jet);% legend('f124fill');
    hold on
   plot(min(.10,f1t )+0.30);
  colormap(jet); legend('f1 f124fill f1t');
    hold on

%--------------------------
  load svff.matJuly24_21C
% load svff.matFeb1_21C
  f24ttfill = fill2D2to1(f24tt);
  f4fill = fill2D4to1(f4);
  f124fill = f12t + fac*f24ttfill + fac*fac*f4fill;

  subplot(2,2,3); 
   plot(min(.10,f1 )+.00);
  colormap(jet);% legend('f1');
    hold on
   plot(min(.10,f124fill )+0.15);
  colormap(jet);% legend('f124fill');
    hold on
   plot(min(.10,f1t )+0.30);
  colormap(jet); legend('f1 f124fill f1t');
    hold on

%--------------------------
  load svff.matJuly24_21D
% load svff.matFeb1_21D
  f24ttfill = fill2D2to1(f24tt);
  f4fill = fill2D4to1(f4);
  f124fill = f12t + fac*f24ttfill + fac*fac*f4fill;

  subplot(2,2,4); 
   plot(min(.10,f1 )+.00);
  colormap(jet);% legend('f1');
    hold on
   plot(min(.10,f124fill )+0.15);
  colormap(jet);% legend('f124fill');
    hold on
   plot(min(.10,f1t )+0.30);
  colormap(jet); legend('f1 f124fill f1t');
    hold on
  
%--------------------------
    hold off
