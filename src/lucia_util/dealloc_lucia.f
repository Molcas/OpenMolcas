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
      SUBROUTINE DEALLOC_LUCIA()
      use stdalloc, only: mma_deallocate
      use GLBBAS, only: INT1, INT1O, PINT1, PINT2, PGINT1, PGINT1A,
     &                  LSM1, LSM2, RHO1, SRHO1, KINH1_NOCCSYM, KINH1,
     &                  CI_VEC, SIGMA_VEC
* Deallocate memory allocated during alloc_lucia

*. Input
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "lucinp.fh"
#include "orbinp.fh"
#include "cstate.fh"
#include "csm.fh"
#include "crun.fh"
#include "cprnt.fh"
#include "cicisp.fh"
*.CSMPRD
#include "csmprd.fh"
#include "cintfo.fh"
*.Output

*.1 : One electron integrals( Complete matrix allocated )
      CALL mma_deallocate(INT1)
*. A copy of the original UNMODIFIED 1-elecs ints
      CALL mma_deallocate(INT1O)
*.2 : Two electron integrals
*. Pointers to symmetry block of integrals
      CALL mma_deallocate(PINT1)
      CALL mma_deallocate(PINT2)
*. Pointers to nonsymmetric one-electron integrals
      DO ISM = 1, NSMOB
*. triangular packed
        CALL mma_deallocate(PGINT1(ISM)%I)
*. no packing
        CALL mma_deallocate(PGINT1A(ISM)%I)
      END DO
*. Symmetry of last index as a function of initial index
      CALL mma_deallocate(LSM1)
      CALL mma_deallocate(LSM2)
*.3 One-body density
      CALL mma_deallocate(RHO1)
*.3.1 : One-body spin density
      CALL mma_deallocate(SRHO1)
*. indices for pair of orbitals symmetry ordered
*. Lower half packed
      CALL mma_deallocate(KINH1)
*. Complete form
      CALL mma_deallocate(KINH1_NOCCSYM)

* arrays allocated at end of lucia
      LBLOCK = MXSOOB
      LBLOCK = MAX(LBLOCK,LCSBLK)
      LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
      IF(PSSIGN.NE.0.0D0) LBLOCK = INT(2.0D0*XISPSM(IREFSM,1))
      Call mma_deallocate(CI_VEC)
      Call mma_deallocate(SIGMA_VEC)
*
      END
