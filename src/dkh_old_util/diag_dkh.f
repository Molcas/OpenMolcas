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
      SUBROUTINE DIAG_DKH(A,N,EIG,EW,SINV,AUX,IC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N*(N+1)/2),AUX(N,N),SINV(N,N),EIG(N,N),EW(N)
#ifdef MOLPRO
CMR
CMR ATTENTION, THE SCRATCH ARRAY TMP IS NOT PROPERLY ALLOCATED
      DIMENSION TMP(6*N)
CMR
#endif
*
      IJ=0
      DO I=1,N
         DO  J=1,I
            IJ=IJ+1
            AUX(I,J)=A(IJ)
            AUX(J,I)=A(IJ)
         END DO
      END DO
      DO 3 K=1,N
        DO 3 J=1,N
          EIG(K,J)=0.D0
          DO 2 L=1,J
    2       EIG(K,J)=EIG(K,J)+AUX(K,L)*SINV(L,J)
    3 CONTINUE
      DO 5 I=1,N
         DO 5 J=1,I
         AUX(I,J)=0.D0
         DO 4 K=1,I
    4      AUX(I,J)=AUX(I,J)+SINV(K,I)*EIG(K,J)
         AUX(J,I)=AUX(I,J)
    5 CONTINUE
*
      TOL=1.D-80
#ifdef MOLPRO
      DO 6 I=1,N
         DO 6 J=1,N
            EIG(J,I)=AUX(J,I)
    6 CONTINUE
      CALL diag2(n,n,ew,eig)
c     call dsyev_('V','L',N,EIG,N,EW,TMP,6*N,INFO)
#else
      CALL JACOB_REL(AUX,EIG,EW,N,TOL,IC)
#endif
*
      RETURN
      END
