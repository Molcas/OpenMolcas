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
      SUBROUTINE MKFOP(FIFA,NGRP,JSTATE_OFF,FOPXMS)

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"

      DIMENSION FIFA(NFIFA)
      DIMENSION FOPXMS(NGRP,NGRP)

* Procedure for computing the matrix FOPXMS, which is the matrix
* of the average Fock operator in the basis of input CASSCF w.f.
* In: The average Fock matrix, active indices only, over the original
* CASSCF MO-orbitals.
* Compute the matrix elements of this Fock operator over the basis
* of CASSCF states that are treated together in XMS fashion.

      CALL QENTER('MKFOP')

      CALL DCOPY_(NGRP**2,[0.0D0],0,FOPXMS,1)

* Loop over bra functions:
      DO J1=1,NGRP
        IBRA=JSTATE_OFF+J1
* Loop over ket functions:
        DO J2=1,J1
          IKET=JSTATE_OFF+J2
* Compute matrix element and put it into FOPXMS:
          CALL FOPAB(FIFA,IBRA,IKET,FOPXMS(J1,J2))
* FOPXMS is symmetric
          FOPXMS(J2,J1) = FOPXMS(J1,J2)
        END DO
      END DO

      CALL QEXIT('MKFOP')
      RETURN
      END
