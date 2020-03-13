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
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE PSGMDIA(ALPHA,BETA,IVEC,JVEC)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

#include "SysDef.fh"

C Compute |JVEC> := BETA*|JVEC> + ALPHA*(H0(diag)-E0)*|IVEC>
C If SHIFT.ne.0.0d0 or SHIFTI.ne.0.0d0, use a modified H0

      DO 100 ICASE=1,13
        DO 100 ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 100
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) GOTO 100
C Remember: NIN values in BDIAG, but must read NAS for correct
C positioning.
          CALL GETMEM('BD','ALLO','REAL',LBD,NAS)
          CALL GETMEM('ID','ALLO','REAL',LID,NIS)
          ID=IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,WORK(LBD),NAS,ID)
          CALL DDAFILE(LUSBT,2,WORK(LID),NIS,ID)

          CALL RHS_ALLO (NIN,NIS,lg_V2)

          IF(BETA.NE.0.0D0) THEN
            CALL RHS_READ (NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
            IF(BETA.NE.1.0D00) THEN
              CALL RHS_SCAL (NIN,NIS,lg_V2,BETA)
            END IF
*         ELSE
*           CALL RHS_SCAL (NIN,NIS,lg_V2,0.0D0)
          END IF

          IF(ALPHA.NE.0.0D0) THEN
            IF(BETA.NE.0.0D0) THEN
              CALL RHS_ALLO (NIN,NIS,lg_V1)
              CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
              CALL RHS_SGMDIA (NIN,NIS,lg_V1,WORK(LBD),WORK(LID))
              CALL RHS_DAXPY(NIN,NIS,ALPHA,lg_V1,lg_V2)
              CALL RHS_FREE (NAS,NIS,lg_V1)
            ELSE
              CALL RHS_READ (NIN,NIS,lg_V2,ICASE,ISYM,IVEC)
              CALL RHS_SGMDIA (NIN,NIS,lg_V2,WORK(LBD),WORK(LID))
              CALL RHS_SCAL (NIN,NIS,lg_V2,ALPHA)
            END IF
          END IF

          CALL RHS_SAVE (NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
          CALL RHS_FREE (NAS,NIS,lg_V2)
          CALL GETMEM('BD','FREE','REAL',LBD,NAS)
          CALL GETMEM('ID','FREE','REAL',LID,NIS)
 100  CONTINUE
      RETURN
      END
