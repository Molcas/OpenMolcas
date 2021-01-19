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
      SUBROUTINE CHO_INIRSDIM()
C
C     Purpose: initialize reduced set dimension.
C
      use ChoArr, only: nDimRS, iSP2F
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRSh, InfRed, IndRed
#include "implicit.fh"
#include "cholesky.fh"

      IF (RSTCHO) THEN
         ILOC = 3
         DO IRS = 1,XNPASS
            CALL CHO_GETRED(INFRED,nnBstRSh(:,:,ILOC),
     &                      IndRed(1,ILOC),INDRSH,iSP2F,
     &                      MAXRED,NSYM,NNSHL,MMBSTRT,IRS,.FALSE.)
            CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,NSYM,NNSHL,ILOC)
            CALL ICOPY(NSYM,NNBSTR(1,ILOC),1,nDimRS(:,iRS),1)
         END DO
         NSET = XNPASS
      ELSE
         CALL ICOPY(NSYM,NNBSTR(1,1),1,nDimRS,1)
         NSET = 1
      END IF

      NUM   = NSYM*(MAXRED - NSET)
      CALL CHO_IZERO(nDimRS(1,NSET+1),NUM)

      END
