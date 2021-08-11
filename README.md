PET Codes
=========

Keith Miller, Aug 6, 2021 </br>
kmiller@math.berkeley.edu

These are the Octave codes used in Linux to produce the figures for my paper

    A Backprojection Kernel (KRNL3D) for Very-Wide-Aperture 3D
    Tomography Applied to PET with Multigrid for Precise Use of
    Time-of_Flight Data

submitted to the journal Physics in Medicine and Biology.

The subdirectories labeled as Fig6, 8, 9, 10, 11, 12 correspond to figures
5, 6, 7, 8, 9, and 10 in the paper. This is because in the final submission
of the paper I deleted the old figures 5 and 7 in order to shorten it.

In each subdirectory `Fig6codesAug_21`, `Fig8codesAug_21`, ... are the codes:
1. `mkTBLK*.m` (make table of values of the kernel K(y1,y2,phi)).
   It calls the subroutines:
   * `H.m` (computes the kernel H(t) for 2D tomography)
   * `vecHH*.m` (computes HH1, HH2, HH4 for blur radius 1,2, or 4)
     It generates `svTBLK1`, or `svTBLK1t`, or `svTBLK124`, etc., containing
     the kernels K1, or K1t, or K124, etc.

2. `randPET*.m` (randomized PET image reconstruction).
  It loads the `svTBLK*` files and calls the subroutines:
  * `onebeam4cube*.m` computes on-the-fly a single (src,dir) with a random
    unit direction vector dir and a source src chosen randomly inside a
    2aaX2bbX2cc "cube". (One could also do an ellipsoid.) That src is
    shifted randomly along the beam to the shifted source "ssrc". With TOF
    info that ssrc becomes "psrc" the "presumed source" from the TOF data.
    The sources are distributed randomly amongst the 4 cubes according to
    the relative frequencies N1, N2, N3, N4 X Nsrc.
  * `vecmap*.m` sets up an orthonormal triplet of vectors with which to
    project the "screen coordinates" y1*, y2* and ssy1, ssy2 of psrc and of
    the reconstruction grid.
  * `lookup*.m` table lookup to get the K(y1,y2,phi) values. The
    reconstructed image matrix is saved in "svff.mat" and then visualized
    by 2D slices using the Octave commands "mesh", "plot", or "contour".

3. view1,1t,or124... can be used to quickly re-view the images f1 or f1t
   or f124 ... in many ways. This after a long time-consuming run.

4. mk2X2newFig6.... this assembles the 4 sub-figures Fig6ABCD for the paper.
   I then annotated, and then used "print -depsc2 Fig6..." to get the .eps
   figures. But WARNING, the .eps versions are HUGE, and needed to be
   compressed - from perhaps 10MB in .eps format to .5Mb in the final .pdf
   format for the paper. (How? my IT expert friend Brandon did it for me.)

NOTE: These codes are extremely NOT-optimized. For example, with multigrid
I should be computing f24tt, f4 on dxx = 2,4 grids. Instead, I compute
them all on the same dxx=1 grid. But then I use the `fill3D2to1.m` and
`fill3D4to1.m` to interpolate, from the coarse-grid values only, onto the
dxx=1 finer grid. Also for example, instead of backprojecting onto the
"backplacement" cylinders, I backproject the truncated kernels K12t, K24tt
with all their zeros onto the entire grid. The success of this project will
require an **efficient** backplacement implementation. But Octave is not the
setting to attempt it. It should probably be done in half-precision on GPUs.
I have included this as item #1 in the list of "needed further work" which
will not be done by me as I return to my retirement (since 2005).

NOTE: I have left dates (such as "mkTBLK124jan29_21.m" on most of these
code names. This is of no use for you, but it helped me to go back to my
handwritten daily notes for info.

NOTE: some labels in the paper, such as LG1t or TOFu, may appear as RG1t
or RTOF in some of the codes.

NOTE: I believe these Octave codes will run seamlessly in Matlab.

NOTE: Feel free to contact me by email if you are serious about running
or adapting these codes and you have questions.
