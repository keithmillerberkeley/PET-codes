


      %  file testfill4to1.m

      close all
      clear all
%  ttt = [1 3 2 2 0 1 0 2 3];
%  ttt = [1 3 2 2 0 1 0 2 3 1 2];
%  ttt = [1 3 2 2 0 1 0 2 3 1 2 0];
   ttt = [1 3 2 2 0 1 0 2 3 1 2 0 2];
   
   sss = fill2D4to1(ttt)

   figure(2);plot(ttt)
   figure(3);plot(sss)

