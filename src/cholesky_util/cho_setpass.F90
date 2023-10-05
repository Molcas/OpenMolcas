!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine CHO_SETPASS(DIAG,DIASH,ISYSH,IRED,CONV,NPOTSH)
!
! Purpose: Check convergence and, if not converged, set up
!          integral pass.

use Cholesky, only: DiaMax, DiaMin, nnShl, nSym, Span, ThrCom
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: Diag(*)
real(kind=wp), intent(_OUT_) :: DIASH(*)
integer(kind=iwp), intent(_OUT_) :: ISYSH(*)
integer(kind=iwp), intent(in) :: IRED
logical(kind=iwp), intent(out) :: CONV
integer(kind=iwp), intent(out) :: NPOTSH
integer(kind=iwp) :: ISHLAB, ISYM
real(kind=wp) :: DGMAX

! Initialize the potential number of shell pairs that can
! contribute.
! -------------------------------------------------------

NPOTSH = 0

! Find max. abs. diagonal in each symmetry and the global max.
! ------------------------------------------------------------

DGMAX = -1.0e15_wp
call CHO_MAXABSDIAG(DIAG,IRED,DGMAX)

! If not converged, set next integral pass.
! -----------------------------------------

CONV = DGMAX < THRCOM
if (.not. CONV) then
  call CHO_SETMAXSHL(DIAG,DIASH,ISYSH,IRED)
  do ISYM=1,NSYM
    DIAMIN(ISYM) = max(DIAMAX(ISYM)*SPAN,THRCOM)
  end do
  do ISHLAB=1,NNSHL
    if (DIASH(ISHLAB) > THRCOM) then
      NPOTSH = NPOTSH+1
    else
      DIASH(ISHLAB) = Zero
    end if
  end do
end if

end subroutine CHO_SETPASS
