************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2017, Roland Lindh                                     *
************************************************************************
* Version of Oct 21

      SUBROUTINE COPDIA(A,VEC,NDIM,IPACK)
*
* Copy diagonal of matrix A into vector VEC
*
*   IPACK = 0 : Full matrix
*   IPACK = 1 : Lower triangular matrix
*
      IMPLICIT REAL*8 ( A-H,O-Z)
      DIMENSION A(*),VEC(*)
      INTEGER   ip_CPDIA
*
#include "WrkSpc.fh"
*VVP
*Workaround for dummy aliasing
      CALL GETMEM('CPDIA','ALLO','REAL',ip_CPDIA,NDIM)
      IF(IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          Work(ip_CPDIA+I-1) = A((I-1)*NDIM+I)
  100   CONTINUE
      ELSE
        DO 200 I = 1, NDIM
          Work(ip_CPDIA+I-1) = A(I*(I+1)/2)
  200   CONTINUE
      END IF
*
      CALL DCOPY_(NDIM,Work(ip_CPDIA),1,VEC,1)
      CALL GETMEM('CPDIA','FREE','REAL',ip_CPDIA,NDIM)

      RETURN
      END
      SUBROUTINE SQRTMT(A,NDIM,ITASK,ASQRT,AMSQRT,SCR)
*
* Calculate square root of positive definite symmetric matrix A
* if(ITASK .EQ. 2 ) Inverted square root matrix is also calculated
* In case of singularities in A A -1/2 is defined to have the same
* singularity
      IMPLICIT REAL*8( A-H,O-Z)
*
      DIMENSION A(NDIM,NDIM)
      DIMENSION ASQRT(NDIM,NDIM),AMSQRT(NDIM,NDIM)
      DIMENSION SCR(*)
* Length of SCR should at least be 2 * NDIM ** 2 + NDIM*(NDIM+1)/2
      KLFREE = 1
*
      KLASYM = KLFREE
      KLAVAL = KLASYM
      KLFREE = KLASYM + NDIM*(NDIM+1)/2
*
      KLAVEC = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
      NTEST = 0
*
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
      CALL TRIPAK(A,SCR(KLASYM),1,NDIM,NDIM)
      Call DCopy_(NDIM**2,[0.0D0],0,SCR(KLAVEC),1)
      Call DCopy_(NDIM,[1.0D0],0,SCR(KLAVEC),1+NDIM)
      Call NIDiag(SCR(KLASYM),SCR(KLAVEC),NDIM,NDIM,0)
      Call JACORD(SCR(KLASYM),SCR(KLAVEC),NDIM,NDIM)
      CALL COPDIA(SCR(KLASYM),SCR(KLAVAL),NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of matrix : '
        CALL WRTMAT(SCR(KLAVAL),NDIM,1,NDIM,1)
      END IF
*. Check for negative eigenvalues
      DO I = 1, NDIM
       IF(Abs(SCR(KLAVAL-1+I)).LT.1.0D-14) SCR(KLAVAL)=0.0D0
       IF(SCR(KLAVAL-1+I).LT.0.0D0) THEN
*         WRITE(6,*) ' SQRTMT : Negative eigenvalue ', SCR(KLAVAL-1+I)
*         WRITE(6,*) ' SQRTMT : I will STOP '
*         STOP       ' SQRTMT : Negative eigenvalue '
        CALL SYSABENDMSG('lucia_util/sqrtmt','Internal error',
     &                   'Negative eigenvalue')
       END IF
      END DO
*
      DO 100 I = 1,NDIM
        SCR(KLAVAL-1+I) = SQRT(SCR(KLAVAL-1+I))
  100 CONTINUE
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(ASQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
*
      IF(ITASK .EQ. 2 ) THEN
        DO 200 I = 1,NDIM
          IF(SCR(KLAVAL-1+I) .GT. 1.0D-13 ) then
            SCR(KLAVAL-1+I) = 1.0D0/SCR(KLAVAL-1+I)
          ELSE
            SCR(KLAVAL-1+I) = SCR(KLAVAL-1+I)
          END IF
  200   CONTINUE
        CALL XDIAXT(AMSQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
      END IF
*
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Info from SQRTMT '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' Input matrix to SQRTMT '
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' Square root of matrix '
        CALL WRTMAT(ASQRT,NDIM,NDIM,NDIM,NDIM)
        IF(ITASK .EQ. 2 ) THEN
          WRITE(6,*) ' Inverse square root of matrix '
          CALL WRTMAT(AMSQRT,NDIM,NDIM,NDIM,NDIM)
        END IF
      END IF
*
      RETURN
      END

      SUBROUTINE XDIAXT(XDX,X,DIA,NDIM,SCR)
*
* Obtain XDX = X * DIA * X(Transposed)
* where DIA is an diagonal matrix stored as a vector
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XDX(NDIM,NDIM)
      DIMENSION X(NDIM,NDIM),DIA(NDIM)
      DIMENSION SCR(NDIM,NDIM)
*
* DIA * X(transposed)
      DO 100 I=1,NDIM
        CALL COPVEC(X(1,I),SCR(1,I),NDIM)
        CALL SCALVE(SCR(1,I),DIA(I),NDIM)
  100 CONTINUE
* X * DIA * X(Transposed)
      CALL MATML4(XDX,X,SCR,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,2)
*
      RETURN
      END
