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
      SUBROUTINE TRACI_RPT2(ISTART,NDIM,XMAT,LSYM,NCI,CI)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XMAT(NDIM,NDIM)
#include "pt2_guga.fh"
#include "WrkSpc.fh"

      IF (NDIM.LE.0) GOTO 999
      NDIM2=NDIM**2

      CALL GETMEM('XSAV','ALLO','REAL',LXSAV,NDIM2)
      CALL DCOPY_(NDIM2,XMAT,1,WORK(LXSAV),1)
      CALL GETMEM('TVEC','ALLO','REAL',LTVEC,NDIM)
      CALL GETMEM('SGM','ALLO','REAL',LSGM,NCI)
      CALL DCOPY_(NCI,0.0D0,0,WORK(LSGM),1)

      DO 100 J=1,NDIM
        FACT=1.0D0/XMAT(J,J)
        DO I=1,NDIM
          WORK(LTVEC-1+I)=-FACT*XMAT(I,J)
          XMAT(I,J)=0.0D0
        END DO
        WORK(LTVEC-1+J)=FACT
        XMAT(J,J)=1.0D00
C Array T now contains a factor of XMAT of the form
C (e(1),..e(k-1),T,..,e(n)), where e(i) is the standard
C unit column vector, with elements Kronecker(l,i).
C Apply its inverse to XMAT.
        DO M=J+1,NDIM
          XJM=XMAT(J,M)
          DO I=1,NDIM
            XMAT(I,M)=XMAT(I,M)+WORK(LTVEC-1+I)*XJM
          END DO
          XMAT(J,M)=WORK(LTVEC-1+J)*XJM
        END DO
C Transform CI array:
C CI:=( 1 + Sum(U(I)E(IJ)) + (1/2)Sum(U(I)U(M)E(IJ,MJ)) ) CI,
C where U(I) = T(I)-Kronecker(I,J).
        JORB=ISTART-1+J
        LJ=LEVEL(JORB)
        CALL DYAX(NCI,(1.5D0-0.5D0*WORK(LTVEC-1+J)),CI,1,WORK(LSGM),1)
        DO I=1,NDIM
          IORB=ISTART-1+I
          LI=LEVEL(IORB)
          SCL=0.5D0*WORK(LTVEC-1+I)
          IF(I.EQ.J) SCL=SCL-0.5D00
          CALL SIGMA1_CP2(LI,LJ,SCL,LSYM,CI,WORK(LSGM),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        END DO
        DO I=1,NDIM
          IORB=ISTART-1+I
          LI=LEVEL(IORB)
          SCL=WORK(LTVEC-1+I)
          IF(I.EQ.J) SCL=SCL-1.0D00
          CALL SIGMA1_CP2(LI,LJ,SCL,LSYM,WORK(LSGM),CI,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        END DO

 100  CONTINUE


      CALL GETMEM('SGM','FREE','REAL',LSGM,NCI)
      CALL GETMEM('TVEC','FREE','REAL',LTVEC,NDIM)
      CALL DCOPY_(NDIM2,WORK(LXSAV),1,XMAT,1)
      CALL GETMEM('XSAV','FREE','REAL',LXSAV,NDIM2)

 999  CONTINUE

      RETURN
      END
