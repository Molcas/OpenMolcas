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
      SUBROUTINE ALLOC_LUCIA()
      use stdalloc, only: mma_allocate
      use GLBBAS, only: INT1, INT1O, PINT1, PINT2, PGINT1, PGINT1A,
     &                  LSM1, LSM2, RHO1, SRHO1, KINH1_NOCCSYM, KINH1
*
* Dimensions and
* Allocation of static memory
*
* =====
* Input
* =====
*
* KFREE : Pointer to first element of free space
* Information in /LUCINP/,/ORBINP/,/CSYM/
*
* ======
* Output
* ======
* KFREE : First array of free space after allocation of
*         static memory
* /GLBBAS/,/CDIM/
*
*
* =======
* Version
* =======
*
* Modified Jan 1997
*           Fall 97 (PGINT1 added )
*           Spring 99

*. Input
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "lucinp.fh"
#include "orbinp.fh"
#include "cstate.fh"
#include "csm.fh"
#include "crun.fh"
#include "cprnt.fh"
*.CSMPRD
#include "csmprd.fh"
#include "cintfo.fh"
*.Output

*.1 : One electron integrals( Complete matrix allocated )
      CALL mma_allocate(INT1,NTOOB ** 2,Label='INT1')
*. A copy of the original UNMODIFIED 1-elecs ints
      CALL mma_allocate(INT1O,NTOOB ** 2,Label='Int1O')
*. Zero to avoid problems with elements that will not
*. be initialized
      ZERO = 0.0D0
      INT1(:)=ZERO
      INT1O(:)=ZERO
*.2 : Two electron integrals
*. Pointers to symmetry block of integrals
      CALL mma_allocate(PINT1,NBINT1,Label='PINT1')
      CALL mma_allocate(PINT2,NBINT2,Label='PINT2')
*. Pointers to nonsymmetric one-electron integrals
      DO ISM = 1, NSMOB
*. triangular packed
        CALL mma_allocate(PGINT1(ISM)%I,NSMOB,Label='PGINT1(ISM)%I')
*. no packing
        CALL mma_allocate(PGINT1A(ISM)%I,NSMOB,Label='PGINT1A(ISM)%I')
      END DO
*. Symmetry of last index as a function of initial index
      CALL mma_allocate(LSM1,NBINT1,Label='LSM1')
      CALL mma_allocate(LSM2,NBINT2,Label='LSM2')
*.3 One-body density
      CALL mma_allocate(RHO1,NTOOB ** 2,Label='RHO1')
*.3.1 : One-body spin density
      CALL mma_allocate(SRHO1,NTOOB **2,Label='SRHO1')
*. indices for pair of orbitals symmetry ordered
*. Lower half packed
      CALL mma_allocate(KINH1,NTOOB*NTOOB,Label='KINTH1')
*. Complete form
      CALL mma_allocate(KINH1_NOCCSYM,NTOOB*NTOOB,
     &                  Label='KINTH1_NOCCSYM')
*
      RETURN
      END
