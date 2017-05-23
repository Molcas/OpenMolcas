************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine MINV(ARRAY,ARRINV,ISING,DET,NDIM)
      Implicit Real*8 (a-h,o-z)
      Real*8 ARRAY(NDIM,NDIM), ARRINV(NDIM,NDIM)
#include "WrkSpc.fh"
*
      Call Allocate_Work(ipA,NDIM**2)
      Call Allocate_Work(ipB,NDIM**2)
      Call Allocate_Work(ipBUF,NDIM)
      Call Allocate_iWork(IPIV,NDIM)
      Call Allocate_iWork(JPIV,NDIM)
*
      Call MINV_(ARRAY,ARRINV,ISING,DET,NDIM,Work(ipA),
     &           Work(ipBUF),Work(ipB),iWork(IPIV),
     &           iWork(JPIV))
*
      Call Free_iWork(JPIV)
      Call Free_iWork(IPIV)
      Call Free_Work(ipBUF)
      Call Free_Work(ipB)
      Call Free_Work(ipA)
*
      Return
      End
      Subroutine MINV_(ARRAY,ARRINV,ISING,DET,NDIM,A,BUF,
     &                 B,IPIV,JPIV)
*     Subroutine Dool(NDIM,MDIM,N,M,A,B,DET,IPIV,JPIV,BUF)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
*
*     SOLVES THE MATRIX EQUATION AX=B BY DOOLITTLE'S METHOD
*     ACTUAL DIMENSIONS ARE N*N AND N*M
*     ALLOCATED DIMENSIONS ARE NDIM*NDIM AND NDIM*MDIM
*     A AND B ARE DESTROYED, AND X IS RETURNED AS MATRIX B
*                                       (MALMQUIST 82-11-12)
*                                        (UPDATE 83-09-28)
*
*                      --- global variables ---
*
      Real*8 ARRAY(NDIM,NDIM), ARRINV(NDIM,NDIM)
      Real*8 A(NDIM,NDIM), BUF(NDIM), B(NDIM,NDIM)
      Integer   IPIV(NDIM), JPIV(NDIM)
*
*     EQUATION IS SOLVED BY FACTORIZING A=L*R IN SAME SPACE AS A.
*     PIVOTING IS ACHIEVED BY INDIRECT INDEXING.
*     FIRST PIVOTING INDICES ARE ASSIGNED START VALUES.
*
*     Call qEnter('MInv')
*
*     set N and M to NDIM since this subroutine is modified to only
*     deal with square matrices.
*
      N = NDIM
      M = NDIM
*     MDIM = NDIM
*
*     Move ARRAY to allocation A and set the B matrix equal to the
*     unity matrix
*
      Do 999 I = 1, NDIM
         Do 998 J = 1, NDIM
            A(I,J) = ARRAY(I,J)
            B(I,J) = Zero
998      Continue
         B(I,I) = One
999   Continue
*
*     Lets go!!
*
      ip = -1
      jp = -1
      Do 1 I=1,N
      IPIV(I)=I
 1    JPIV(I)=I
      DET=One
      Do 5 I=1,N
*
*     NOW FIND BETTER PIVOT ELEMENT:
*
      AMAX=-One
      Do 2 K=I,N
      Do 2 L=I,N
      AM=ABS(A(IPIV(K),JPIV(L)))
      IF(AMAX.GT.AM) Go To 2
      AMAX=AM
      IP=K
      JP=L
 2    Continue
      IF(IP.eq.I) Go To 3
      DET=-DET
      IDUM=IPIV(I)
      IPIV(I)=IPIV(IP)
      IPIV(IP)=IDUM
 3    IF(JP.eq.I) Go To 4
      DET=-DET
      IDUM=JPIV(I)
      JPIV(I)=JPIV(JP)
      JPIV(JP)=IDUM
 4    IP=IPIV(I)
      JP=JPIV(I)
      DIAG=A(IP,JP)
      BUF(I)=DIAG
      DET=DET*DIAG
      Do 5 K=I+1,N
      KP=IPIV(K)
      C=A(KP,JP)
      IF(DIAG.NE.Zero) C=C/DIAG
      A(KP,JP)=C
      Do 5 L=I+1,N
      LP=JPIV(L)
 5    A(KP,LP)=A(KP,LP)-C*A(IP,LP)
*
*     FIRST RESUBSTITUTION STEP:
*
      Do 7 J=1,M
      Do 7 I=2,N
      IP=IPIV(I)
      SUM=B(IP,J)
      Do 6 K=1,I-1
 6    SUM=SUM-A(IP,JPIV(K))*B(IPIV(K),J)
 7    B(IP,J)=SUM
*
*     SECOND RESUBSTITUTION STEP:
*
      Do 9 J=1,M
      Do 9 I=N,1,-1
      IP=IPIV(I)
      SUM=B(IP,J)
      Do 8 K=I+1,N
 8    SUM=SUM-A(IP,JPIV(K))*B(IPIV(K),J)
      IF(BUF(I).NE.Zero) SUM=SUM/BUF(I)
 9    B(IP,J)=SUM
*
*     REORGANIZATION PART:
*
      Do 12 J=1,M
      Do 10 I=1,N
10    BUF(I)=B(IPIV(I),J)
      Do 11 I=1,N
11    B(JPIV(I),J)=BUF(I)
12    Continue
*
*     Move the result to location ARRINV
*
      Do 888 I = 1, NDIM
         Do 889 J = 1, NDIM
            ARRINV (I,J) = B(I,J)
889      Continue
888   Continue
*
*     Call qExit('MInv')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(ISING)
      End
