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

implicit none
#include "chcc1.fh"
real*8 W2(1:nv,1:nv*(nv+1)/2)
integer a, be, ga, bega, bad
real*8 s

bad = 0

bega = 0
do be=1,nv
  do ga=1,be
    bega = bega+1

    do a=1,nv

      s = Q4(a,ga,a,be)/2
      if (abs(W2(a,bega)-s) > 1.0d-10) bad = bad+1
      W2(a,bega) = s

    end do

  end do
end do

write(6,*) ' W2+ chck ',bad

return

end subroutine Chck_W2p
