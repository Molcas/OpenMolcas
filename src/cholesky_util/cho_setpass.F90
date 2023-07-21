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

implicit real*8(a-h,o-z)
real*8 Diag(*), DIASH(*)
integer ISYSH(*)
logical CONV
#include "cholesky.fh"

! Initialize the potential number of shell pairs that can
! contribute.
! -------------------------------------------------------

NPOTSH = 0

! Find max. abs. diagonal in each symmetry and the global max.
! ------------------------------------------------------------

DGMAX = -1.0d15
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
      DIASH(ISHLAB) = 0.0d0
    end if
  end do
end if

end subroutine CHO_SETPASS
