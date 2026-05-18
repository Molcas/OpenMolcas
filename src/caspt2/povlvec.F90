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
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! New VECTOR UTILITIES, written by Steven Vancoillie, May 2011
! A set of subroutines that can transform RHS arrays using the parallel
! aware subroutines
!***********************************************************************

      SUBROUTINE POVLVEC (IVEC,JVEC,OVLAPS)
      use constants, only: Zero
      use caspt2_module, only: NINDEP, NISUP, MxCase, nCASES, CPUOVL,   &
     &                         TIOOVL, nSym
      use definitions, only: wp, iwp
      IMPLICIT None

      real(kind=wp), intent(out) :: OVLAPS(0:8,0:MXCASE)
      integer(kind=iwp), intent(in) :: iVec, jVec

      real(kind=wp) CPU, CPU0, CPU1
      real(kind=wp) TIO, TIO0, TIO1
      real(kind=wp) OVLTOT, OVL, OVLSUM
      integer(kind=iwp) iCase, iSym, lg_v1, lg_v2, NIN, NIS
      real(kind=wp), External :: RHS_DDOT

! Compute overlaps of vectors nr IVEC and JVEC in SR format!, for each
! individual case and symmetry block, in OVLAPS(ISYM,ICASE), summed over
! symmetry in OVLAPS(0,ICASE), summed over case in OVLAPS(ISYM,0), total
! sum in OVLAPS(0,0).
      CALL TIMING(CPU0,CPU,TIO0,TIO)

      OVLTOT=Zero
      OVLAPS(:,0)=Zero

      DO ICASE=1,NCASES
        OVLSUM=Zero
        DO ISYM=1,NSYM
          OVL=Zero
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NIN*NIS.NE.0) THEN
            CALL RHS_ALLO (NIN,NIS,lg_V1)
            CALL RHS_READ (NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
            IF (IVEC.NE.JVEC) THEN
              CALL RHS_ALLO (NIN,NIS,lg_V2)
              CALL RHS_READ (NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
            ELSE
              lg_V2=lg_V1
            END IF
            OVL=RHS_DDOT(NIN,NIS,lg_V1,lg_V2)
            CALL RHS_FREE (lg_V1)
            IF (IVEC.NE.JVEC) THEN
              CALL RHS_FREE (lg_V2)
            END IF
          END IF
          OVLAPS(ISYM,ICASE)=OVL
          OVLAPS(ISYM,0)=OVLAPS(ISYM,0)+OVL
          OVLSUM=OVLSUM+OVL
        END DO
        OVLAPS(0,ICASE)=OVLSUM
        OVLTOT=OVLTOT+OVLSUM
      END DO
      OVLAPS(0,0)=OVLTOT

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUOVL=CPUOVL+(CPU1-CPU0)
      TIOOVL=TIOOVL+(TIO1-TIO0)

      END SUBROUTINE POVLVEC
