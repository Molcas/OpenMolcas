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
      SUBROUTINE CHO_RSTOF(IRS2F,N,LRDIM,IRED)
C
C     Purpose: set up mapping between reduced set and SO indices
C              (i.e., full storage).
C
C     IRS2F(1,irs) = alpha (SO index, not symmmetry reduced)
C     IRS2F(2,irs) = beta  (SO index, not symmmetry reduced)
C
#include "implicit.fh"
      INTEGER IRS2F(N,LRDIM)
#include "cholesky.fh"
#include "choorb.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_RSTOF')

      INTEGER  CHO_RS2F
      EXTERNAL CHO_RS2F
      INTEGER  CHO_F2SP
      EXTERNAL CHO_F2SP

      MULD2H(I,J)=IEOR(I-1,J-1)+1
      ITRI(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
      ISOSHL(I)=IWORK(ip_iSOShl-1+I)
      ISHLSO(I)=IWORK(ip_iShlSO-1+I)
      NBSTSH(I)=IWORK(ip_NBSTSH-1+I)

      IF (N .LT. 2) THEN
         CALL CHO_QUIT('Dimension error [1] in '//SECNAM,104)
      END IF
      IF (LRDIM .NE. MMBSTRT) THEN
         CALL CHO_QUIT('Dimension error [2] in '//SECNAM,104)
      END IF
      CALL CHO_IZERO(IRS2F,N*MMBSTRT)

      DO ISYMA = 1,NSYM
         IF (NBAS(ISYMA) .GT. 0) THEN
            DO ISYMB = 1,ISYMA-1
               ISYMAB = MULD2H(ISYMA,ISYMB)
               DO KB = 1,NBAS(ISYMB)
                  IB = IBAS(ISYMB) + KB
                  LB = ISHLSO(IB)
                  ISHLB = ISOSHL(IB)
                  DO KA = 1,NBAS(ISYMA)
                     IA = IBAS(ISYMA) + KA
                     LA = ISHLSO(IA)
                     ISHLA  = ISOSHL(IA)
                     IF (ISHLA .LT. ISHLB) THEN
                        LAB = NBSTSH(ISHLB)*(LA - 1) + LB
                     ELSE IF (ISHLA .EQ. ISHLB) THEN
                        LAB = ITRI(LA,LB)
                     ELSE
                        LAB = NBSTSH(ISHLA)*(LB - 1) + LA
                     END IF
                     ISHLAB = CHO_F2SP(ITRI(ISHLA,ISHLB))
                     IF (ISHLAB .GT. 0) THEN
                        IRS = CHO_RS2F(LAB,ISHLAB,ISYMAB,IRED)
                        IF (IRS .GT. 0) THEN
                           IRS2F(1,IRS) = IA
                           IRS2F(2,IRS) = IB
                        END IF
                     END IF
                  END DO
               END DO
            END DO
            ISYMB  = ISYMA
            ISYMAB = 1
            DO KA = 1,NBAS(ISYMA)
               IA = IBAS(ISYMA) + KA
               LA = ISHLSO(IA)
               ISHLA  = ISOSHL(IA)
               DO KB = 1,KA
                  IB = IBAS(ISYMB) + KB
                  LB = ISHLSO(IB)
                  ISHLB = ISOSHL(IB)
                  IF (ISHLA .LT. ISHLB) THEN
                     LAB = NBSTSH(ISHLB)*(LA - 1) + LB
                  ELSE IF (ISHLA .EQ. ISHLB) THEN
                     LAB = ITRI(LA,LB)
                  ELSE
                     LAB = NBSTSH(ISHLA)*(LB - 1) + LA
                  END IF
                  ISHLAB = CHO_F2SP(ITRI(ISHLA,ISHLB))
                  IF (ISHLAB .GT. 0) THEN
                     IRS = CHO_RS2F(LAB,ISHLAB,ISYMAB,IRED)
                     IF (IRS .GT. 0) THEN
                        IRS2F(1,IRS) = IA
                        IRS2F(2,IRS) = IB
                     END IF
                  END IF
               END DO
            END DO
         END IF
      END DO

      END
