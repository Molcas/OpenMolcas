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

subroutine DerHeffX(IVEC,JVEC,NASHT,NTG3,OVL,DTG1,DTG2,DTG3)
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

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: GA_Arrays
use caspt2_module, only: NASUP, NINDEP, NISUP, NSYM
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC, NASHT, NTG3
real(kind=wp), intent(out) :: OVL
real(kind=wp), intent(inout) :: DTG1(NASHT,NASHT), DTG2(NASHT,NASHT,NASHT,NASHT), DTG3(NTG3)
integer(kind=iwp) :: ICASE, iHi1, iHi2, iLo1, iLo2, ISYM, jHi1, jHi2, jLo1, jLo2, lg_V1, lg_V2, MV1, MV2, NAS, NIN, NIS, nvlen
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

! The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

! Sketch of procedure:
!  Loop over every (case/symmetry)-block.
!     If (No such vector block) Skip to end of loop
!     Allocate two places for this block, VEC1 and VEC2
!     Read VEC1 as IVEC component from file.
!     Read VEC2 as JVEC component from file.
!     Loop nest, computing
!        HEL := HEL + VEC1*GOM*VEC2
!     End of loop nest
!     Deallocate VEC1 and VEC2
!  End of loop.

do ICASE=1,13
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIN = NINDEP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)

    if (NAS*NIS == 0) cycle
    if (NIN == 0) cycle

    call RHS_ALLO(NAS,NIS,lg_V1)
    call RHS_ALLO(NAS,NIS,lg_V2)
    call RHS_READ(NAS,NIS,lg_V1,ICASE,ISYM,IVEC)
    call RHS_READ(NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
    call RHS_ACCESS(NAS,NIS,lg_V1,iLo1,iHi1,jLo1,jHi1,MV1)
    call RHS_ACCESS(NAS,NIS,lg_V2,iLo2,iHi2,jLo2,jHi2,MV2)

    if ((iLo1 /= iLo2) .or. (iHi1 /= iHi2) .or. (jLo1 /= jLo2) .or. (jHi1 /= jHi2)) then
      write(u6,'(1X,A)') 'HCOUP: Error: block mismatch, abort...'
      call ABEND()
    end if

    nvlen = (iHi1-jLo1+1)*(jHi1-jLo1+1)

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call DerHEffX_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,jLo1,jHi1,DBL_MB(MV1),DBL_MB(MV2),OVL,DTG1,DTG2,DTG3)
    else
#   endif
      call DerHEffX_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,jLo1,jHi1,GA_Arrays(MV1)%A,GA_Arrays(MV2)%A,OVL,DTG1,DTG2,DTG3)
#   ifdef _MOLCAS_MPP_
    end if
#   endif

    call RHS_RELEASE(lg_V1,iLo1,iHi1,jLo1,jHi1)
    call RHS_RELEASE(lg_V2,iLo2,iHi2,jLo2,jHi2)
    call RHS_FREE(lg_V1)
    call RHS_FREE(lg_V2)
  end do
end do

end subroutine DerHeffX
