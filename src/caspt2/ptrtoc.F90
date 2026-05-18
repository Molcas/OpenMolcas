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

subroutine PTRTOC(ITYPE,IVEC,JVEC)
! Transform RHS vectors from SR format to C format.
! ITYPE=0 uses only T matrix, ITYPE=1 uses S*T matrix

use caspt2_module, only: CPUVEC, NASUP, nCases, nInDep, NISUP, nSym, TIOVEC
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iType, iVec, jVec
real(kind=wp) :: CPU, CPU0, CPU1, TIO, TIO0, TIO1
integer(kind=iwp) :: iCase, iSym, lg_v1, lg_v2, NAS, NIN, NIS

call TIMING(CPU0,CPU,TIO0,TIO)

do ICASE=1,NCASES
  if ((IVEC == JVEC) .and. (ICASE == 12)) cycle
  if ((IVEC == JVEC) .and. (ICASE == 13)) cycle
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    if (NAS == 0) cycle
    NIN = NINDEP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NIS == 0) cycle
    call RHS_ALLO(NAS,NIS,lg_V2)
    if ((ICASE /= 12) .and. (ICASE /= 13)) then
      if (NIN > 0) then
        call RHS_ALLO(NIN,NIS,lg_V1)
        call RHS_READ(NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
        call RHS_SR2C(ITYPE,0,NAS,NIS,NIN,lg_V1,lg_V2,ICASE,ISYM)
        call RHS_FREE(lg_V1)
      else
        call RHS_SCAL(NAS,NIS,lg_V2,Zero)
      end if
    else
      call RHS_READ(NAS,NIS,lg_V2,ICASE,ISYM,IVEC)
    end if

    call RHS_SAVE(NAS,NIS,lg_V2,ICASE,ISYM,JVEC)

    call RHS_FREE(lg_V2)

  end do
end do

call TIMING(CPU1,CPU,TIO1,TIO)
CPUVEC = CPUVEC+(CPU1-CPU0)
TIOVEC = TIOVEC+(TIO1-TIO0)

end subroutine PTRTOC
