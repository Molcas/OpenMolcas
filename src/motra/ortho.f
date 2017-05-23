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
      SUBROUTINE ORTHO_MOTRA(nSym,nBas,nDel,Ovlp,CMO)
*
*     Objective: Orthonormalize input vectors
*               (Gram-Schmidt orthogonaliztion)
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "trafo_motra.fh"
#include "WrkSpc.fh"
*
      Real*8 Ovlp(*), CMO(*)
      Integer nBas(nSym), nDel(nSym)
*
      Call qEnter('ORTHO')
*
*     Allocate work space
*
      CALL GETMEM('SCR1','ALLO','REAL',LW1,N2MAX)
      CALL GETMEM('SCR2','ALLO','REAL',LW2,N2MAX)
      CALL GETMEM('SCR3','ALLO','REAL',LW3,N2MAX)
*
*     Loop over symmetries and orthogonalize
*
      IJ=1
      IM=1
      DO ISYM=1,NSYM
         NORBI=NBAS(ISYM)-NDEL(ISYM)
         IF(NORBI.GT.0) THEN
           CALL SQUARE(OVLP(IJ),WORK(LW3),1,NBAS(ISYM),NBAS(ISYM))
C           CALL MXMA(WORK(LW3),1,NBAS(ISYM),
C     *               CMO(IM),1,NBAS(ISYM),
C     *               WORK(LW2),1,NBAS(ISYM),
C     *               NBAS(ISYM),NBAS(ISYM),NORBI)
           CALL DGEMM_('N','N',NBAS(ISYM),NORBI,NBAS(ISYM),
     *                  1.0d0,WORK(LW3),NBAS(ISYM),CMO(IM),
     *                  NBAS(ISYM),0.0d0,WORK(LW2),NBAS(ISYM))
C           CALL MXMA(CMO(IM),NBAS(ISYM),1,
C     *               WORK(LW2),1,NBAS(ISYM),
C     *               WORK(LW1),1,NORBI,
C     *               NORBI,NBAS(ISYM),NORBI)
           CALL DGEMM_('T','N',NORBI,NORBI,NBAS(ISYM),
     *                  1.0d0,CMO(IM),NBAS(ISYM),WORK(LW2),
     *                  NBAS(ISYM),0.0d0,WORK(LW1),NORBI)
           CALL ORTHOX_MOTRA(WORK(LW1),CMO(IM),NORBI,NBAS(ISYM))
        END IF
           IJ=IJ+NBAS(ISYM)*(NBAS(ISYM)+1)/2
           IM=IM+NBAS(ISYM)*NBAS(ISYM)
      END DO
*
*     Deallocate work space and exit
*
      CALL GETMEM('SCR3','FREE','REAL',LW3,N2MAX)
      CALL GETMEM('SCR2','FREE','REAL',LW2,N2MAX)
      CALL GETMEM('SCR1','FREE','REAL',LW1,N2MAX)
*
      Call qExit('ORTHO')
*
      RETURN
      END
