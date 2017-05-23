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
      SUBROUTINE FZERO(B,N)
      INTEGER    N
      REAL*8     B(N)

      CALL DCOPY_(N,0.0D0,0,B,1)

      RETURN
      END

      SUBROUTINE IZERO(B,N)
      INTEGER    N
      INTEGER    B(N)

      CALL ICOPY(N,0,0,B,1)

      RETURN
      END

      SUBROUTINE DFILL (N,SA,SX,INCX)
      INTEGER    N,INCX
      REAL*8     SX(N),SA

      CALL DCOPY_(N,SA,0,SX,INCX)

      RETURN
      END


      SUBROUTINE FMOVE(IA,IB,N)
C      INTEGER    N
C      REAL*8     IA(N),IB(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION   IA(*),IB(*)
#include "SysDef.fh"

C      CALL DCOPY_(N*RtoI,IA,1,IB,1)
      DO I=1,N*RtoI
         IB(I)=IA(I)
      END DO

      RETURN
      END

      SUBROUTINE ADDVEC (A,B,C,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N),B(N),C(N)

      DO I=1,N
         A(I)=B(I)+C(I)
      END DO

      RETURN
      END

      SUBROUTINE SUBVEC (A,B,C,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N),B(N),C(N)

      DO I=1,N
         A(I)=B(I)-C(I)
      END DO

      RETURN
      END


      SUBROUTINE SCATTER(N,A,IND,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(N),IND(N)

      DO I=1,N
         A(IND(I))=B(I)
      END DO

      RETURN
      END
