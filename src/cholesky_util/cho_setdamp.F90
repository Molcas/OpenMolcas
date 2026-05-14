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

subroutine CHO_SETDAMP()
!
! Purpose: set screening damping, unless user-defined.

use Cholesky, only: Damp, ThrCom
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: I

do I=1,2
  if (DAMP(I) < Zero) then
    if (THRCOM > 9.99e-3_wp) then ! >= 1.0e-2
      DAMP(I) = 1.0e7_wp
    else if (THRCOM > 9.99e-4_wp) then ! >= 1.0e-3
      DAMP(I) = 1.0e6_wp
    else if (THRCOM > 9.99e-5_wp) then ! >= 1.0e-4
      DAMP(I) = 1.0e5_wp
    else if (THRCOM > 9.99e-6_wp) then ! >= 1.0e-5
      DAMP(I) = 1.0e4_wp
    else if (THRCOM > 9.99e-7_wp) then ! >= 1.0e-6
      DAMP(I) = 1.0e3_wp
    else if (THRCOM > 9.99e-8_wp) then ! >= 1.0e-7
      DAMP(I) = 1.0e2_wp
    else if (THRCOM > 9.99e-9_wp) then ! >= 1.0e-8
      DAMP(I) = 1.0e1_wp
    else
      DAMP(I) = One
    end if
  end if
end do

end subroutine CHO_SETDAMP
