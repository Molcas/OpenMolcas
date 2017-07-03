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
      SUBROUTINE SETSXCI

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "output_ras.fh"
      PARAMETER(ROUTINE='SETSXCI ')

      Common /IDSXCI/ IDXCI(mxAct),IDXSX(mxAct)

      DIMENSION IOFF_GSSH(mxgas)
C
C ---------------------------------------------------------
C --  SET INDEX VECTORS FOR CI/SX INTEGRAL ORDERING
C ---------------------------------------------------------
C
      Call QENTER(ROUTINE)

      NGSSHT=0
      DO IGAS=1,NGAS
        IOFF_GSSH(IGAS)=NGSSHT
        NGSSHT=NGSSHT+SUM(NGSSH(IGAS,1:NSYM))
      END DO
      ISTOT=0
      DO ISYM=1,NSYM
        DO IGAS=1,NGAS
          DO IGSSH=1,NGSSH(IGAS,ISYM)
            IOFF_GSSH(IGAS)=IOFF_GSSH(IGAS)+1
            ISTOT=ISTOT+1
            IDXCI(ISTOT)=IOFF_GSSH(IGAS)
          END DO
        END DO
      END DO

      DO I=1,ISTOT
        IDXSX(IDXCI(I))=I
      END DO

      IF (IPRGLB.GE.DEBUG)THEN
        WRITE(6,'(1X,A,1X,12I5)')
     &   'REORDERING VECTOR FOR CI',(IDXCI(I),I=1,ISTOT)
        WRITE(6,'(1X,A,1X,12I5)')
     &   'REORDERING VECTOR FOR SX',(IDXSX(I),I=1,ISTOT)
      ENDIF
      RETURN
      END
