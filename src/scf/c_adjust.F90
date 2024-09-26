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

subroutine C_Adjust(CInter,n,CThr)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: CInter(n), CThr
integer(kind=iwp) :: i
real(kind=wp) :: Fact

if (CInter(n) < CThr) then
  Fact = (One-CThr)/(One-CInter(n))
  do i=1,n-1
    CInter(i) = Fact*CInter(i)
  end do
  CInter(n) = CThr
end if

return

end subroutine C_Adjust
