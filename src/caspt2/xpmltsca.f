************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
#include "xrhs.fh"
      SUBROUTINE PMLTSCA(KOD,IMLTOP,LST1,LST2,
     &                   X,NXI,NXA,F,NFI,NFA,
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
      DIMENSION X(NXI,NXA),F(NFI,NFA)
      DIMENSION LST1(4,NLST1), LST2(4,NLST2)

#ifdef _MOLCAS_MPP_
C SVC: Determine the index ranges of the local chunks of lg_Y.
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
          IF (KOD.EQ.23 .OR. KOD.EQ.24) THEN
            CALL MLTSCA_DH(IMLTOP,LST1,LST2,
     &                   X,NXI,NXA,F,NFI,NFA,
     &                   DBL_MB(mY),NAS2,jYLo,jYHi)
          ELSE
            WRITE(6,*) 'PMLTSCA: not supposed to be here'
            CALL AbEnd()
          END IF
          CALL GA_Release_Update (lg_Y,iYLo,iYHi,jYLo,jYHi)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        IF (KOD.EQ.23 .OR. KOD.EQ.24) THEN
          CALL MLTSCA_DH(IMLTOP,LST1,LST2,
     &                   X,NXI,NXA,F,NFI,NFA,
     &                   GA_Arrays(lg_Y)%Array,NAS2,1,NIS2)
        ELSE
          WRITE(6,*) 'PMLTSCA: not supposed to be here'
          CALL AbEnd()
        END IF
#ifdef _MOLCAS_MPP_
      END IF
#endif
      END SUBROUTINE PMLTSCA
