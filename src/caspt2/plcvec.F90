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

subroutine PLCVEC(ALPHA,BETA,IVEC,JVEC)

use constants, only: Zero, One
use caspt2_module, only: nCases, nSym, nInDep, NISUP, CPULCS, TIOLCS
use definitions, only: iwp, wp

implicit none

real(kind=wp), intent(in) :: Alpha, Beta
integer(kind=iwp), intent(in) :: iVec, jVec
real(kind=wp) CPU, CPU0, CPU1
real(kind=wp) TIO, TIO0, TIO1
integer(kind=iwp) iCase, iSym, lg_v1, lg_v2
integer(kind=iwp) NIN, NIS

! |JVEC> := BETA*|JVEC> + ALPHA*|IVEC>, IVEC and JVEC in SR format!

if ((BETA == One) .and. (ALPHA == Zero)) return

call TIMING(CPU0,CPU,TIO0,TIO)

do ICASE=1,NCASES
  do ISYM=1,NSYM
    NIN = NINDEP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NIN*NIS == 0) cycle

    call RHS_ALLO(NIN,NIS,lg_V2)
    if ((BETA /= Zero) .and. (ALPHA /= Zero)) call RHS_ALLO(NIN,NIS,lg_V1)

    if ((BETA == Zero) .and. (ALPHA == Zero)) then
      call RHS_SCAL(NIN,NIS,lg_V2,Zero)
    else if (BETA /= Zero) then
      call RHS_READ(NIN,NIS,lg_V2,ICASE,ISYM,JVEC)
      if (ALPHA /= Zero) then
        call RHS_SCAL(NIN,NIS,lg_V2,BETA)
        call RHS_READ(NIN,NIS,lg_V1,ICASE,ISYM,IVEC)
        call RHS_DAXPY(NIN,NIS,ALPHA,lg_V1,lg_V2)
      else
        call RHS_SCAL(NIN,NIS,lg_V2,BETA)
      end if
    else if (ALPHA /= Zero) then
      call RHS_READ(NIN,NIS,lg_V2,ICASE,ISYM,IVEC)
      call RHS_SCAL(NIN,NIS,lg_V2,ALPHA)
    end if

    call RHS_SAVE(NIN,NIS,lg_V2,ICASE,ISYM,JVEC)

    call RHS_FREE(lg_V2)
    if ((BETA /= Zero) .and. (ALPHA /= Zero)) call RHS_FREE(lg_V1)
  end do
end do

call TIMING(CPU1,CPU,TIO1,TIO)
CPULCS = CPULCS+(CPU1-CPU0)
TIOLCS = TIOLCS+(TIO1-TIO0)

end subroutine PLCVEC
