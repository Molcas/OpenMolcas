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

subroutine Chck_W2p(W2)
! test W2+

use Index_Functions, only: nTri_Elem
use chcc_global, only: nv, Q4
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: W2(nv,nTri_Elem(nv))
integer(kind=iwp) :: a, bad, be, bega, ga
real(kind=wp) :: s

bad = 0

bega = 0
do be=1,nv
  do ga=1,be
    bega = bega+1

    do a=1,nv

      s = Q4(a,ga,a,be)*Half
      if (abs(W2(a,bega)-s) > 1.0e-10_wp) bad = bad+1
      W2(a,bega) = s

    end do

  end do
end do

write(u6,*) ' W2+ chck ',bad

return

end subroutine Chck_W2p
