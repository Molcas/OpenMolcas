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
      SUBROUTINE DEALLOC_LUCIA
* Deallocate memory allocated during alloc_lucia

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
#include "cicisp.fh"
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
      CALL GETMEM('INT1  ','FREE','REAL',KINT1,NTOOB ** 2)
*. A copy of the original UNMODIFIED 1-elecs ints
      CALL GETMEM('INT1O ','FREE','REAL',KINT1O,NTOOB ** 2)
*.1.1 : Inactive fock matrix
      CALL GETMEM('FI    ','FREE','REAL',KFI  ,NTOOB ** 2)
      CALL GETMEM('FIO   ','FREE','REAL',KFIO ,NTOOB ** 2)
*.1.2 Inactive Fock matrx in zero order space
      CALL GETMEM('FIZ    ','FREE','REAL',KFIZ,NTOOB **2)
*.2 : Two electron integrals
*. Pointers to symmetry block of integrals
      CALL GETMEM('PINT1 ','FREE','INTE',KPINT1,NBINT1)
      CALL GETMEM('PINT2 ','FREE','INTE',KPINT2,NBINT2)
*. Pointers to nonsymmetric one-electron integrals
      DO ISM = 1, NSMOB
*. triangular packed
        CALL GETMEM('PGINT1','FREE','INTE',KPGINT1(ISM),NSMOB)
*. no packing
        CALL GETMEM('PGIN1A','FREE','INTE',KPGINT1A(ISM),NSMOB)
      END DO
*. Symmetry of last index as a function of initial index
      CALL GETMEM('LSM1   ','FREE','INTE',KLSM1,NBINT1)
      CALL GETMEM('LSM2   ','FREE','INTE',KLSM2,NBINT2)
*.3 One-body density
      CALL GETMEM('RHO1  ','FREE','REAL',KRHO1,NTOOB ** 2)
*.3.1 : One-body spin density
      CALL GETMEM('SRHO1 ','FREE','REAL',KSRHO1,NTOOB **2)
*. indices for pair of orbitals symmetry ordered
*. Lower half packed
      CALL GETMEM('KINH1  ','FREE','INTE',KINH1,NTOOB*NTOOB)
*. Complete form
      CALL GETMEM('KINH1  ','FREE','INTE',KINH1_NOCCSYM,NTOOB*NTOOB)

* arrays allocated at end of lucia
      LBLOCK = MXSOOB
      LBLOCK = MAX(LBLOCK,LCSBLK)
      LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
      IF(PSSIGN.NE.0.0D0) LBLOCK = INT(2.0D0*XISPSM(IREFSM,1))
      CALL GETMEM('VEC1  ','FREE','REAL',KCI_POINTER,LBLOCK)
      CALL GETMEM('VEC2  ','FREE','REAL',KSIGMA_POINTER,LBLOCK)
*
      RETURN
      END
