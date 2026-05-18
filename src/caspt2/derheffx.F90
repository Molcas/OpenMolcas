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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
      SUBROUTINE DerHeffX(IVEC,JVEC,NASHT,NTG3,OVL,DTG1,DTG2,DTG3)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: NSYM, NASUP, NISUP, NINDEP
      use definitions, only: wp, iwp, u6

      implicit none
! Compute the coupling Hamiltonian element defined as
!     HEL = < ROOT1 | H * OMEGA | ROOT2 >
! assuming that IVEC contains a contravariant representation of
! H|ROOT1>, JVEC contains a contravariant representation of
! OMEGA|ROOT2>, and OVL, TG1, TG2, TG3 contain the overlap (normally
! expected to be 0 or 1) and active transition density matrices of ROOT1
! and ROOT2. See also subroutine TSVEC for explanations.

! SVC (March 2014): modification of original code to handle distributed
! RHS arrays. There is now a main HCOUP subroutine that loops over cases
! and irreps and gets access to the process-specific block of the RHS.
! The coupling for that block is computed by the subroutine HCOUP_BLK.

      integer(kind=iwp), intent(in) :: IVEC, JVEC, NASHT, NTG3
      real(kind=wp), intent(out) :: OVL
      real(kind=wp), intent(inout) :: DTG1(NASHT,NASHT),                &
     &  DTG2(NASHT,NASHT,NASHT,NASHT), DTG3(NTG3)

      integer(kind=iwp) :: ICASE, ISYM, NAS, NIN, NIS, lg_V1, lg_V2,    &
     &  iLo1, iHi1, jLo1, jHi1, MV1, iLo2, iHi2, jLo2, jHi2, MV2
      integer(kind=iwp) :: nvlen
! The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif


! Sketch of procedure:
!  Loop over every (case/symmetry)-block.
!           If (No such vector block) Skip to end of loop
!           Allocate two places for this block, VEC1 and VEC2
!           Read VEC1 as IVEC component from file.
!           Read VEC2 as JVEC component from file.
!           Loop nest, computing
!              HEL := HEL + VEC1*GOM*VEC2
!           End of loop nest
!           Deallocate VEC1 and VEC2
!  End of loop.

      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)

          IF(NAS*NIS == 0) cycle
          IF(NIN == 0) cycle

          CALL RHS_ALLO (NAS,NIS,lg_V1)
          CALL RHS_ALLO (NAS,NIS,lg_V2)
          CALL RHS_READ (NAS,NIS,lg_V1,ICASE,ISYM,IVEC)
          CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
          CALL RHS_ACCESS(NAS,NIS,lg_V1,iLo1,iHi1,jLo1,jHi1,MV1)
          CALL RHS_ACCESS(NAS,NIS,lg_V2,iLo2,iHi2,jLo2,jHi2,MV2)

          IF ((iLo1 /= iLo2) .OR. (iHi1 /= iHi2) .OR.                   &
     &        (jLo1 /= jLo2) .OR. (jHi1 /= jHi2)) THEN
            WRITE(u6,'(1X,A)') 'HCOUP: Error: block mismatch, abort...'
            CALL ABEND()
          END IF

          nvlen = (iHi1-jLo1+1)*(jHi1-jLo1+1)

#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            CALL DerHEffX_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,jLo1,jHi1,&
     &                      DBL_MB(MV1),DBL_MB(MV2),OVL,                &
     &                      DTG1,DTG2,DTG3)
          ELSE
#endif
            CALL DerHEffX_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,jLo1,jHi1,&
     &                        GA_Arrays(MV1)%A,                         &
     &                        GA_Arrays(MV2)%A,OVL,                     &
     &                        DTG1,DTG2,DTG3)
#ifdef _MOLCAS_MPP_
          END IF
#endif

          CALL RHS_RELEASE (lg_V1,iLo1,iHi1,jLo1,jHi1)
          CALL RHS_RELEASE (lg_V2,iLo2,iHi2,jLo2,jHi2)
          CALL RHS_FREE (lg_V1)
          CALL RHS_FREE (lg_V2)
        END DO
      END DO

      end subroutine DerHeffX
