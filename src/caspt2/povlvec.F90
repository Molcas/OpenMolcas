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

subroutine POVLVEC(IVEC,JVEC,OVLAPS)
! Compute overlaps of vectors nr IVEC and JVEC in SR format!, for each
! individual case and symmetry block, in OVLAPS(ISYM,ICASE), summed over
! symmetry in OVLAPS(0,ICASE), summed over case in OVLAPS(ISYM,0), total
! sum in OVLAPS(0,0).

use caspt2_module, only: CPUOVL, MxCase, nCASES, NINDEP, NISUP, nSym, TIOOVL
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iVec, jVec
real(kind=wp), intent(out) :: OVLAPS(0:8,0:MXCASE)
integer(kind=iwp) :: iCase, iSym, lg_v1, lg_v2, NIN, NIS
real(kind=wp) :: CPU, CPU0, CPU1, OVL, OVLSUM, OVLTOT, TIO, TIO0, TIO1
real(kind=wp), external :: RHS_DDOT

call TIMING(CPU0,CPU,TIO0,TIO)

OVLTOT = Zero
OVLAPS(:,0) = Zero

do ICASE=1,NCASES
  OVLSUM = Zero
  do ISYM=1,NSYM
    OVL = Zero
    NIN = NINDEP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NIN*NIS /= 0) then
      call RHS_ALLO(NIN,NIS,lg_V1)
      call RHS_READ(NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
      if (IVEC /= JVEC) then
        call RHS_ALLO(NIN,NIS,lg_V2)
        call RHS_READ(NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
      else
        lg_V2 = lg_V1
      end if
      OVL = RHS_DDOT(NIN,NIS,lg_V1,lg_V2)
      call RHS_FREE(lg_V1)
      if (IVEC /= JVEC) call RHS_FREE(lg_V2)
    end if
    OVLAPS(ISYM,ICASE) = OVL
    OVLAPS(ISYM,0) = OVLAPS(ISYM,0)+OVL
    OVLSUM = OVLSUM+OVL
  end do
  OVLAPS(0,ICASE) = OVLSUM
  OVLTOT = OVLTOT+OVLSUM
end do
OVLAPS(0,0) = OVLTOT

call TIMING(CPU1,CPU,TIO1,TIO)
CPUOVL = CPUOVL+(CPU1-CPU0)
TIOOVL = TIOOVL+(TIO1-TIO0)

end subroutine POVLVEC
