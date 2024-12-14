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
* Copyright (C) 1991,1997, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE INTIM()
      use GLBBAS, only: INT1, INT1O, PINT1, PINT2, LSM1, LSM2, KINH1,
     &                  KINH1_NOCCSYM
      use lucia_data, only: ECORE_HEX,ECORE_ORIG,ECORE
      use lucia_data, only: NOINT
*
* Interface to external integrals
*
* If NOINT .ne. 0, only pointers are constructed
* Jeppe Olsen, Winter of 1991
*
* Version : Fall 97
*
      IMPLICIT None
#include "mxpdim.fh"
#include "lucinp.fh"
#include "csm.fh"
#include "orbinp.fh"
#include "oper.fh"
*
*. : Pointers for symmetry blocks of integrals
*
      CALL INTPNT(PINT1,LSM1,PINT2,LSM2)
*
*. Pointer for orbital indices for symmetry blocked matrices
      CALL ORBINH1(KINH1,KINH1_NOCCSYM,NTOOBS,NTOOB,NSMOB)
*
*. Change one-electron integrals to inactive fock matrix
      IF(NOINT.EQ.0) THEN
        INT1O(:)=INT1(:)
        ECORE_HEX = 0.0D0
      END IF
      ECORE_ORIG = ECORE
      ECORE = ECORE + ECORE_HEX

      END SUBROUTINE INTIM
