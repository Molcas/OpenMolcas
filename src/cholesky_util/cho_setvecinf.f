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
      SUBROUTINE CHO_SETVECINF(INFVEC,MVEC,M2,MSYM,IVEC,ISYM,IAB,IPASS,
     &                         ILOC)
C
C     Purpose: set info for vector IVEC of sym. ISYM.
C
#include "implicit.fh"
      INTEGER INFVEC(MVEC,M2,MSYM)
#include "cholesky.fh"

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_SETVECINF')

#if defined (_DEBUGPRINT_)
      CALL QENTER('_SETVECINF')
#endif

      IF (IVEC .GT. MAXVEC) THEN
         WRITE(LUPRI,*) SECNAM,': too many Cholesky vectors!'
         WRITE(LUPRI,*) SECNAM,': symmetry: ',ISYM
         WRITE(LUPRI,*) SECNAM,': max. allowed is ',MAXVEC
         WRITE(LUPRI,*) SECNAM,': please increase max. ',
     &                  'allowed'
         CALL CHO_QUIT('Too many Cholesky vectors in '
     &                 //SECNAM,104)
      ELSE IF (IVEC .EQ. MAXVEC) THEN ! no set next addr.
         INFVEC(IVEC,1,ISYM) = IAB    ! diag. index red. set 1
         INFVEC(IVEC,2,ISYM) = IPASS  ! global red. set
      ELSE
         INFVEC(IVEC,1,ISYM)   = IAB   ! diag. index red. set 1
         INFVEC(IVEC,2,ISYM)   = IPASS ! global red. set
         INFVEC(IVEC+1,4,ISYM) = INFVEC(IVEC,4,ISYM)
     &                         + NNBSTR(ISYM,ILOC) ! next addr.
      END IF

#if defined (_DEBUGPRINT_)
      CALL QEXIT('_SETVECINF')
#endif

      END
