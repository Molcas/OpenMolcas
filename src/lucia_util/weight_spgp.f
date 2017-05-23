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
      SUBROUTINE WEIGHT_SPGP(      Z, NORBTP, NELFTP,NORBFTP,   ISCR,
     &                         NTEST)
*
* construct vertex weights for given supergroup
*
* Reverse lexical ordering is used
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*. Input
      INTEGER NELFTP(NORBTP),NORBFTP(NORBTP)
*. Ouput
      INTEGER Z(*)
*. Scratch length : 2 * NORB + (NEL+1)*(NORB+1)
       INTEGER ISCR(*)
*
       NORB = IELSUM(NORBFTP,NORBTP)
       NEL  = IELSUM(NELFTP,NORBTP)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Subroutine WEIGHT_SPGP in action '
        WRITE(6,*) ' ================================='
        WRITE(6,*) 'NELFTP '
        CALL IWRTMA(NELFTP,1,NORBTP,1,NORBTP)
      END IF
*
      KLFREE = 1
      KLMAX = KLFREE
      KLFREE = KLFREE + NORB
*
      KLMIN = KLFREE
      KLFREE = KLFREE + NORB
*
      KW = KLFREE
      KLFREE = KW + (NEL+1)*(NORB+1)
*.Max and min arrays for strings
      CALL MXMNOC_SPGP(ISCR(KLMIN),ISCR(KLMAX),NORBTP,NORBFTP,NELFTP,
     &                         NTEST)
*. Arc weights
      CALL GRAPW(  ISCR(KW),         Z,ISCR(KLMIN),ISCR(KLMAX),    NORB,
     &                  NEL,     NTEST)
*
      RETURN
      END
