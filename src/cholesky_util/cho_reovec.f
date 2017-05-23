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
      SUBROUTINE CHO_REOVEC(IRS2F,N,LRDIM,WRK,LWRK)
C
C     Purpose: reorder Cholesky vectors on disk to full storage.
C
#include "implicit.fh"
      INTEGER   IRS2F(N,LRDIM)
      DIMENSION WRK(LWRK)
#include "cholesky.fh"
#include "choorb.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_REOVEC')

      INTEGER  CHO_ISAO
      EXTERNAL CHO_ISAO

      ITRI(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

      CALL QENTER('_REOVEC')

C     Set up mapping from rs1 to full storage.
C     ----------------------------------------

      IF (N .LT. 3) THEN
         CALL CHO_QUIT('Dimension error [1] in '//SECNAM,104)
      END IF
      IF (LRDIM .NE. MMBSTRT) THEN
         CALL CHO_QUIT('Dimension error [2] in '//SECNAM,104)
      END IF
      CALL CHO_RSTOF(IRS2F,N,MMBSTRT,1)
      DO IRS1 = 1,NNBSTRT(1)
         IA = IRS2F(1,IRS1)
         IB = IRS2F(2,IRS1)
         ISYMA = CHO_ISAO(IA)
         ISYMB = CHO_ISAO(IB)
         JA = IA - IBAS(ISYMA)
         JB = IB - IBAS(ISYMB)
         IF (ISYMA .EQ. ISYMB) THEN
            JAB = ITRI(JA,JB)
         ELSE
           JAB = NBAS(ISYMA)*(JB - 1) + JA
         END IF
         IRS2F(1,IRS1) = ISYMA
         IRS2F(2,IRS1) = ISYMB
         IRS2F(3,IRS1) = JAB
      END DO

C     Set up index arrays and open files for storing full vectors.
C     ------------------------------------------------------------

      CALL CHO_REOINI()

C     Reorder vectors on disk.
C     ------------------------

      CALL CHO_REOVC1(IRS2F,N,LRDIM,WRK,LWRK)

      CALL QEXIT('_REOVEC')

      END
