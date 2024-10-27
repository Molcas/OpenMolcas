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
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ALPHA, BETA
      INTEGER IVEC, JVEC

#include "caspt2.fh"
      INTEGER ICASE, ISYM, NIN, NAS, NIS, JD
      REAL*8, ALLOCATABLE:: BD(:), ID(:)

C Compute |JVEC> := BETA*|JVEC> + ALPHA*(H0(diag)-E0)*|IVEC>
C If real_shift.ne.0.0d0 or imag_shift.ne.0.0d0, use a modified H0

      DO 100 ICASE=1,13
        DO 101 ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 101
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) GOTO 101
C Remember: NIN values in BDIAG, but must read NAS for correct
C positioning.
          CALL mma_allocate(BD,NAS,LABEL='BD')
          CALL mma_allocate(ID,NIS,LABEL='ID')
          JD=IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,BD,NAS,JD)
          CALL DDAFILE(LUSBT,2,ID,NIS,JD)

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
              CALL RHS_SGMDIA (NIN,NIS,lg_V1,BD,ID)
              CALL RHS_DAXPY(NIN,NIS,ALPHA,lg_V1,lg_V2)
              CALL RHS_FREE (lg_V1)
            ELSE
              CALL RHS_READ (NIN,NIS,lg_V2,ICASE,ISYM,IVEC)
              CALL RHS_SGMDIA (NIN,NIS,lg_V2,BD,ID)
              CALL RHS_SCAL (NIN,NIS,lg_V2,ALPHA)
            END IF
          END IF

          CALL RHS_SAVE (NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
          CALL RHS_FREE (lg_V2)
          CALL mma_deallocate(BD)
          CALL mma_deallocate(ID)
 101    CONTINUE
 100  CONTINUE
      RETURN
      END
