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
      SUBROUTINE WEIGHT_mclr(Z,NEL,NORB1,NORB2,NORB3,
     &                  MNRS1,MXRS1,MNRS3,MXRS3, ISCR)
*
* construct vertex weights
*
* Reverse lexical ordering is used for restricted space
*
      IMPLICIT REAL*8 ( A-H,O-Z)
      INTEGER Z(*), ISCR(*)
      NTEST=0
*
      NORB = NORB1 + NORB2 + NORB3
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
      CALL RSMXMN(ISCR(KLMAX),ISCR(KLMIN),NORB1,NORB2,NORB3,
     &            NEL,MNRS1,MXRS1,MNRS3,MXRS3,NTEST)
*. Arc weights
      CALL GRAPW(ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NEL,NTEST)
*
      RETURN
      END
