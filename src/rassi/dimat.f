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
* Copyright (C) 2000, Per Ake Malmqvist                                *
************************************************************************
*****************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST 2000-06-30
*  SUBROUTINE DIMAT
*  CONSTRUCT A DENSITY MATRIX DINAO FOR THE INACTIVE ORBITALS.
*  IT IS RETURNED IN SYMMETRY-BLOCKED SQUARED FORMAT.
*****************************************************************
      SUBROUTINE DIMAT(CMO1,CMO2,DINAO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CMO1(NCMO),CMO2(NCMO),DINAO(NBSQ)
#include "rassi.fh"
#include "symmul.fh"
      CALL DCOPY_(NBSQ,[0.0D0],0,DINAO,1)
      ISTC=1
      ISTD=1
      DO ISYM=1,NSYM
       NI=NISH(ISYM)
       NO=NOSH(ISYM)
       NB=NBASF(ISYM)
       IF(NI.NE.0) THEN
        CALL  DGEMM_('N','T',NB,NB,NI,1.0D0,
     &               CMO1(ISTC),NB,
     &               CMO2(ISTC),NB,0.0D0,
     &               DINAO(ISTD),NB)
       END IF
       ISTC=ISTC+NO*NB
       ISTD=ISTD+NB**2
      END DO
      CALL DSCAL_(NBSQ,2.0D00,DINAO,1)
      RETURN
      END
