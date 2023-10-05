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

subroutine CHO_CHKINTO(XINT,DIAG,ISYM,NERR,TOL,REPORT)
!
! Purpose: check diagonals in qualified integral columns against
!          original diagonal (read in here).

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: XINT(*), TOL
real(kind=wp), intent(_OUT_) :: DIAG(*)
integer(kind=iwp), intent(in) :: ISYM
integer(kind=iwp), intent(out) :: NERR
logical(kind=iwp), intent(in) :: REPORT
integer(kind=iwp) :: IOPT

IOPT = 2
call CHO_IODIAG(DIAG,IOPT)
call CHO_P_CHKINT(XINT,DIAG,ISYM,NERR,TOL,REPORT)

end subroutine CHO_CHKINTO
