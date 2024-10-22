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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
#include "xrhs.fh"
C
C-----------------------------------------------------------------------
C
      !! RHS_SGMDIA
      SUBROUTINE CASPT2_ResD(Mode,NIN,NIS,lg_W1,lg_W2,DIN,DIS)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use EQSOLV
      use fake_GA, only: GA_Arrays
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
      REAL*8 DIN(*),DIS(*)

C Apply the resolvent of the diagonal part of H0 to an RHS array

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
C-SVC: get the local vertical stripes of the lg_W vector
      !! not yet, probably just repeat for lg_W2
        CALL GA_Distribution (lg_W1,myRank,iLo,iHi,jLo,jHi)
        IF (iLo.NE.0.AND.jLo.NE.0) THEN
          NROW=iHi-iLo+1
          NCOL=jHi-jLo+1
          CALL GA_Access (lg_W1,iLo,iHi,jLo,jHi,mW,LDW)
          CALL CASPT2_ResD2(MODE,NROW,NCOL,DBL_MB(mW),DBL_MB(mW),
     &                      LDW,DIN(iLo),DIS(jLo))
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        CALL CASPT2_ResD2(MODE,NIN,NIS,GA_Arrays(lg_W1)%Array,
     &                                 GA_Arrays(lg_W2)%Array,
     *                                 NIN,DIN,DIS)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE CASPT2_ResD
C
C-----------------------------------------------------------------------
C
      !! RESDIA
