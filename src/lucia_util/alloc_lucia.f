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
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
#include "cintfo.fh"
#include "rasscf_lucia.fh"
*.Output
#include "glbbas.fh"

*.1 : One electron integrals( Complete matrix allocated )
      CALL GETMEM('INT1  ','ALLO','REAL',KINT1,NTOOB ** 2)
*. A copy of the original UNMODIFIED 1-elecs ints
      CALL GETMEM('INT1O ','ALLO','REAL',KINT1O,NTOOB ** 2)
      kint1_pointer = KINT1
      kint1o_pointer = KINT1O
*. Zero to avoid problems with elements that will not
*. be initialized
      ZERO = 0.0D0
      CALL SETVEC(WORK(KINT1),ZERO,NTOOB**2)
      CALL SETVEC(WORK(KINT1O),ZERO,NTOOB**2)
*.1.1 : Inactive fock matrix
      CALL GETMEM('FI    ','ALLO','REAL',KFI  ,NTOOB ** 2)
      CALL GETMEM('FIO   ','ALLO','REAL',KFIO ,NTOOB ** 2)
*.1.2 Inactive Fock matrx in zero order space
      CALL GETMEM('FIZ    ','ALLO','REAL',KFIZ,NTOOB **2)
*.2 : Two electron integrals
*. Pointers to symmetry block of integrals
      CALL GETMEM('PINT1 ','ALLO','INTE',KPINT1,NBINT1)
      CALL GETMEM('PINT2 ','ALLO','INTE',KPINT2,NBINT2)
*. Pointers to nonsymmetric one-electron integrals
      DO ISM = 1, NSMOB
*. triangular packed
        CALL GETMEM('PGINT1','ALLO','INTE',KPGINT1(ISM),NSMOB)
*. no packing
        CALL GETMEM('PGIN1A','ALLO','INTE',KPGINT1A(ISM),NSMOB)
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
