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
* Copyright (C) 2018, Jesper Norell                                    *
************************************************************************
      SUBROUTINE DYSON(IFSBTAB1,
     &    IFSBTAB2,ISSTAB,
     &    DET1,DET2,
     &    IF10,IF01,DYSAMP,DYSCOF)

      IMPLICIT NONE
      INTEGER IFSBTAB1(*),IFSBTAB2(*),ISSTAB(*)
      INTEGER LORBTB
      INTEGER NASORB
      REAL*8 DET1(*),DET2(*),DYSAMP,DYSCOF(*)
      LOGICAL IF10,IF01

#include "symmul.fh"
#include "WrkSpc.fh"

! +++ J.Norell 12/7 - 2018
C Given two CI expansions, using a biorthonormal set of SD's,
C (assuming one state with N and one with N-1 electrons)
C calculate the following quantities:
C (1) The Dyson orbital norm between the two states (DYSAMP)
C (2) The Dyson orbital expansion coefficients in the
C     basis of active biorthonormal orbitals (DYSCOF)
C More functionality should be added here later.

      LORBTB=ISSTAB(3)
C Pick out nr of active orbitals from orbital table:
      NASORB=IWORK(LORBTB+3)

      CALL MKDYSORB(IWORK(LORBTB),ISSTAB,
     &               IFSBTAB1,IFSBTAB2,DET1,DET2,
     &               IF10,IF01,DYSAMP,DYSCOF)

      RETURN
      END
