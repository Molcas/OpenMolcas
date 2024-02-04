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
      SUBROUTINE ALLOC_LUCIA
      use stdalloc, only: mma_allocate
      use GLBBAS
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
*           Fall 97 (KPGINT1 added )
*           Spring 99

*. Input
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "lucinp.fh"
#include "orbinp.fh"
#include "cstate.fh"
#include "csm.fh"
#include "crun.fh"
#include "cprnt.fh"
*.CSMPRD
#include "csmprd.fh"
#include "cintfo.fh"
#include "rasscf_lucia.fh"
*.Output

*.1 : One electron integrals( Complete matrix allocated )
      CALL mma_allocate(INT1,NTOOB ** 2,Label='INT1')
*. A copy of the original UNMODIFIED 1-elecs ints
      CALL mma_allocate(INT1O,NTOOB ** 2,Label='Int1O')
      kint1_pointer = ip_of_Work(INT1(1))
      kint1o_pointer = ip_of_Work(INT1O(1))
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
        CALL GETMEM('PGINT1','ALLO','INTE',KPGINT1(ISM),NSMOB)
*. no packing
        CALL mma_allocate(PGINT1A(ISM)%I,NSMOB,Label='PGINT1A(ISM)%I')
      END DO
*. Symmetry of last index as a function of initial index
      CALL GETMEM('LSM1   ','ALLO','INTE',KLSM1,NBINT1)
      CALL GETMEM('LSM2   ','ALLO','INTE',KLSM2,NBINT2)
*.3 One-body density
      CALL GETMEM('RHO1  ','ALLO','REAL',KRHO1,NTOOB ** 2)
*.3.1 : One-body spin density
      CALL GETMEM('SRHO1 ','ALLO','REAL',KSRHO1,NTOOB **2)
*. indices for pair of orbitals symmetry ordered
*. Lower half packed
      CALL GETMEM('KINH1  ','ALLO','INTE',KINH1,NTOOB*NTOOB)
*. Complete form
      CALL GETMEM('KINH1  ','ALLO','INTE',KINH1_NOCCSYM,NTOOB*NTOOB)
*
      RETURN
      END
