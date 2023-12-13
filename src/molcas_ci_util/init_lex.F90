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

subroutine INIT_LEX(K,LEX)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: K
integer(kind=iwp), intent(out) :: LEX(K)
integer(kind=iwp) :: I

do I=1,K
  LEX(I) = I
end do

end subroutine INIT_LEX
