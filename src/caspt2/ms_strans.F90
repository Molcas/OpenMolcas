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

      SUBROUTINE MS_STRANS(IVEC,JVEC,NASHT,NTG3,OVL,TG1,TG2,TG3,HEL,    &
     &                     SCAL)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: NSYM, NASUP, NISUP, NINDEP, CASES
      use Constants, only: Zero
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
      real(kind=wp), intent(inout) :: OVL, TG1(NASHT,NASHT),            &
     &  TG2(NASHT,NASHT,NASHT,NASHT), TG3(NTG3)
      real(kind=wp), intent(out) :: HEL
      real(kind=wp), intent(in) :: SCAL
! The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

      real(kind=wp) :: HECOMP(14,9), HEBLK, SUMSYM, SUMCASE
      integer(kind=iwp) :: ICASE, ISYM, NAS, NIN, NIS, lg_V1, iLo1,     &
     &  iHi1, jLo1, jHi1, MV1, lg_V2, iLo2, iHi2, jLo2, jHi2, MV2,      &
     &  NHECOMP, i, IC, IS
      integer(kind=iwp) :: nvlen

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

      HEL=Zero
      HECOMP=Zero
      DO ICASE=1,13
!     if (icase /= 12.and.icase /= 13) cycle ! H
!     if (icase /= 10.and.icase /= 11) cycle ! G
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          HEBLK=Zero

          if (NAS*NIS /= 0 .and. NIN /= 0) then
            CALL RHS_ALLO (NAS,NIS,lg_V1)
            CALL RHS_ALLO (NAS,NIS,lg_V2)
            CALL RHS_READ (NAS,NIS,lg_V1,ICASE,ISYM,IVEC)
            CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
            CALL RHS_ACCESS(NAS,NIS,lg_V1,iLo1,iHi1,jLo1,jHi1,MV1)
            CALL RHS_ACCESS(NAS,NIS,lg_V2,iLo2,iHi2,jLo2,jHi2,MV2)

            IF ((iLo1 /= iLo2) .OR. (iHi1 /= iHi2) .OR.                 &
     &          (jLo1 /= jLo2) .OR. (jHi1 /= jHi2)) THEN
              WRITE(u6,'(1X,A)')'HCOUP: Error: block mismatch, abort...'
              CALL ABEND()
            END IF

            nvlen = (iHi1-jLo1+1)*(jHi1-jLo1+1)

#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              CALL MS_STRANS_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,       &
     &                        jLo1,jHi1,DBL_MB(MV1),DBL_MB(MV2),OVL,    &
     &                        TG1,TG2,TG3,SCAL)
            ELSE
#endif
              CALL MS_STRANS_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,       &
     &                        jLo1,jHi1,GA_Arrays(MV1)%A,               &
     &                        GA_Arrays(MV2)%A,OVL,                     &
     &                        TG1,TG2,TG3,SCAL)
#ifdef _MOLCAS_MPP_
            END IF
#endif
            !! Save T*S
            CALL RHS_SAVE (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)

            CALL RHS_RELEASE (lg_V1,iLo1,iHi1,jLo1,jHi1)
            CALL RHS_RELEASE (lg_V2,iLo2,iHi2,jLo2,jHi2)
            CALL RHS_FREE (lg_V1)
            CALL RHS_FREE (lg_V2)
          end if
          HECOMP(ICASE,ISYM)=HEBLK
          HEL=HEL+HEBLK
        END DO
      END DO

! Sum-reduce the per-process contributions
      CALL GADGOP_SCAL(HEL,'+')
      NHECOMP=14*9
      CALL GADGOP(HECOMP,NHECOMP,'+')

      IF(IPRGLB >= DEBUG) THEN
        DO ICASE=1,13
          SUMSYM=Zero
          DO ISYM=1,NSYM
            SUMSYM=SUMSYM+HECOMP(ICASE,ISYM)
          END DO
          HECOMP(ICASE,NSYM+1)=SUMSYM
        END DO

        DO ISYM=1,NSYM+1
          SUMCASE=Zero
          DO ICASE=1,13
            SUMCASE=SUMCASE+HECOMP(ICASE,ISYM)
          END DO
          HECOMP(14,ISYM)=SUMCASE
        END DO

        WRITE(u6,'(20a4)')('----',i=1,20)
        WRITE(u6,*)                                                     &
     &          'HCOUP: The contributions to the Hamiltonian coupling'
        WRITE(u6,*)' elements, by case and by symmetry label.'
        DO IC=1,13
          WRITE(u6,'(1X,A8,9F12.8)')                                    &
     &      CASES(IC),(HECOMP(IC,IS),IS=1,NSYM+1)
        END DO
        CALL XFLUSH(u6)
        WRITE(u6,'(1X,A8,9F12.8)')                                      &
     &    'Summed: ', (HECOMP(14,IS),IS=1,NSYM+1)
        WRITE(u6,*)
      END IF

      end subroutine MS_STRANS
