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
      SUBROUTINE PMLTR1 (KOD,IMLTOP,LST1,
     &                   X,NAS1,NIS1,JXOFF,
     &                   F,NFI,NFJ,
     &                   lg_Y,NAS2,NIS2)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use Sigma_data
      use fake_GA, only: GA_Arrays
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      INTEGER KOD, IMLTOP
      INTEGER LST1(4,NLST1)
      REAL*8 X(*)
      INTEGER NAS1,NIS1,JXOFF, NFI,NFJ, lg_Y, NAS2, NIS2
      REAL*8 F(NFI,NFJ)

#ifdef _MOLCAS_MPP_
C SVC: Determine the index ranges of the local chunks of lg_X and lg_Y.
C The boundaries and leading dimension are stored in a common block for
C access inside the lower-level routines.
C For now, only case H is handled as a distributed array, which is
C always the Y array.
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
*     CALL GA_Distribution (lg_X,myRank,iXLo,iXHi,jXLo,jXHi)
*     IF (iXLo.NE.0.AND.jXLo.NE.0) THEN
*       CALL GA_Access (lg_X,iXLo,iXHi,jXLo,jXHi,mX,LDX)
*     END IF
        CALL GA_Distribution (lg_Y,myRank,iYLo,iYHi,jYLo,jYHi)
        IF (iYLo.NE.0.AND.jYLo.NE.0) THEN
          CALL GA_Access (lg_Y,iYLo,iYHi,jYLo,jYHi,mY,LDY)
          IF (KOD.EQ.17 .OR. KOD.EQ.18) THEN
            CALL MLTR1_EH(IMLTOP,LST1,
     &                  X,NAS1,NIS1,JXOFF,
     &                  F,NFI,NFJ,
     &                  DBL_MB(mY),NAS2,jYLo,jYHi)
          ELSE IF (KOD.EQ.21 .OR. KOD.EQ.22) THEN
            CALL MLTR1_GH(IMLTOP,LST1,
     &                  X,NAS1,NIS1,JXOFF,
     &                  F,NFI,NFJ,
     &                  DBL_MB(mY),NAS2,jYLo,jYHi)
          END IF
          CALL GA_Release_Update (lg_Y,iYLo,iYHi,jYLo,jYHi)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        IF (KOD.EQ.17 .OR. KOD.EQ.18) THEN
          CALL MLTR1_EH(IMLTOP,LST1,
     &                  X,NAS1,NIS1,JXOFF,
     &                  F,NFI,NFJ,
     &                  GA_Arrays(lg_Y)%A,NAS2,1,NIS2)
        ELSE IF (KOD.EQ.21 .OR. KOD.EQ.22) THEN
          CALL MLTR1_GH(IMLTOP,LST1,
     &                  X,NAS1,NIS1,JXOFF,
     &                  F,NFI,NFJ,
     &                  GA_Arrays(lg_Y)%A,NAS2,1,NIS2)
        END IF
#ifdef _MOLCAS_MPP_
      END IF
#endif
      RETURN
      END
