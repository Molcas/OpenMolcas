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

implicit none
#include "cholesky.fh"
integer I

do I=1,2
  if (DAMP(I) < 0.0d0) then
    if (THRCOM > 9.99D-3) then ! >= 1.0D-2
      DAMP(I) = 1.0d7
    else if (THRCOM > 9.99D-4) then ! >= 1.0D-3
      DAMP(I) = 1.0d6
    else if (THRCOM > 9.99D-5) then ! >= 1.0D-4
      DAMP(I) = 1.0d5
    else if (THRCOM > 9.99D-6) then ! >= 1.0D-5
      DAMP(I) = 1.0d4
    else if (THRCOM > 9.99D-7) then ! >= 1.0D-6
      DAMP(I) = 1.0d3
    else if (THRCOM > 9.99D-8) then ! >= 1.0D-7
      DAMP(I) = 1.0d2
    else if (THRCOM > 9.99D-9) then ! >= 1.0D-8
      DAMP(I) = 1.0d1
    else
      DAMP(I) = 1.0d0
    end if
  end if
end do

end subroutine CHO_SETDAMP
