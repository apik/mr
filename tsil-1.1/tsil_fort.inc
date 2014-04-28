c     This is the include file that defines what is necessary for the
c     generic FORTRAN interface with REAL*8 (double) data.

      COMPLEX*16 M,U(4),T(6),S(2),B(2),V(4),Tbar(6)
      COMPLEX*16 Ubold(4,3),Vbold(4,3),Tbold(6,3),Sbold(2,3)

c     This common block matches a C struct defined in the wrapper code
c     fevaluate.c:
      COMMON /results/ M,U,T,S,B,V,Tbar,Ubold,Vbold,Tbold,Sbold

c     This value is used to check for data size compatibility with the 
c     library when tsilfevaluate is called:
      INTEGER rsize
      COMMON /rsize/ rsize
      DATA rsize/8/

c     Variables for indexing arrays:
      INTEGER xz, yu
      PARAMETER(xz=1, yu=2)

      INTEGER vyz, uxv, yzv, xuv, zyv, vxu
      PARAMETER(vyz=1, uxv=2, yzv=3, xuv=4, zyv=5, vxu=6)

      INTEGER zxyv, uyxv, xzuv, yuzv
      PARAMETER(zxyv=1, uyxv=2, xzuv=3, yuzv=4)

