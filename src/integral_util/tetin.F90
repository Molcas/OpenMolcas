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

subroutine tetin(lmax)

use Constants, only: Zero, One
use welcom, only: binom, tetint
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lmax
integer(kind=iwp) :: i, k, l, lm2, m

do k=0,lmax
  lm2 = k/2
  do l=0,lm2
    tetint(k,l) = Zero
    m = k-l*2
    do i=0,l
      tetint(k,l) = tetint(k,l)+binom(l,i)*(-One)**i/real(m+i*2+1,kind=wp)
    end do
  end do
end do

end subroutine tetin
