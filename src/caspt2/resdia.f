!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
!> @brief Apply the resolvent of the diagonal of H0-E0 to a vector
!> @author Per-&Aring;ke Malmqvist, 1994
!>
!> @details
!> ::presdia stands for "parallel resolvent diagonal" and it applies
!> the resolvent of the diagonal of \f$ (H_0 - E_0) \f$ to a coefficient
!> vector stored at position \p IVEC on \p LUSOLV. The result is stored
!> in the vector at position \p JVEC. Furthremore, it computes the
!> "overlap" as \p IVEC squared divided by the resolvent.
!> The parallel part is thaken care of in subroutine calls \p RHS_.
!> Potential shifts and modification of \f$ H_0 \f$ are considered
!> in ::rhs_resdia
!>
!> @param[in]     IVEC   Vector position to which the res is applied
!> @param[out]    JVEC   Vector position where the result is saved
!> @param[in]     OVLAPS Array containing the overlaps
      SUBROUTINE PRESDIA(IVEC,JVEC,OVLAPS)
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

#include "SysDef.fh"

      INTEGER IVEC,JVEC
      REAL*8 OVLAPS(0:8,0:MXCASE)

      INTEGER ICASE,ISYM
      INTEGER NAS,NIS,NIN
      INTEGER LBD,LID,lg_V
      INTEGER ID

      REAL*8 OVL,DOVL,OVLSUM,OVLTOT

C Apply the resolvent of the diagonal part of H0 to a coefficient
C vector in vector nr. IVEC on LUSOLV. Put the results in vector
C nr. JVEC. Also compute overlaps, see OVLVEC for structure.

      OVLTOT=0.0D0
      DO ISYM=1,NSYM
        OVLAPS(ISYM,0)=0.0D0
      END DO
      DO 200 ICASE=1,13
        OVLSUM=0.0D0
        DO 100 ISYM=1,NSYM
          OVL=0.0D0
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) GOTO 90
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
C Remember: NIN values in BDIAG, but must read NAS for correct
C positioning.
          CALL GETMEM('LBD','ALLO','REAL',LBD,NAS)
          CALL GETMEM('LID','ALLO','REAL',LID,NIS)
          ID=IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,WORK(LBD),NAS,ID)
          CALL DDAFILE(LUSBT,2,WORK(LID),NIS,ID)

          CALL RHS_ALLO(NIN,NIS,lg_V)
          CALL RHS_READ(NIN,NIS,lg_V,ICASE,ISYM,IVEC)
          CALL RHS_RESDIA(NIN,NIS,lg_V,WORK(LBD),WORK(LID),DOVL)
          OVL=OVL+DOVL
          CALL RHS_SAVE(NIN,NIS,lg_V,ICASE,ISYM,JVEC)
          CALL RHS_FREE(NIN,NIS,lg_V)

          CALL GETMEM('LBD','FREE','REAL',LBD,NAS)
          CALL GETMEM('LID','FREE','REAL',LID,NIS)
  90      CONTINUE
          OVLAPS(ISYM,0)=OVLAPS(ISYM,0)+OVL
          OVLSUM=OVLSUM+OVL
 100    CONTINUE
        OVLAPS(0,ICASE)=OVLSUM
        OVLTOT=OVLTOT+OVLSUM
 200  CONTINUE
      OVLAPS(0,0)=OVLTOT

      RETURN
      END
