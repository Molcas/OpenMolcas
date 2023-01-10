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

subroutine check_loops(nv,vblock,nla,nlb)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nv, vblock
integer(kind=iwp), intent(out) :: nla, nlb
integer(kind=iwp) :: nga, ngb, ngc, nuga

nuga = nv/vblock
!mp! pridavok
if ((nuga*vblock) < nv) nuga = nuga+1

nla = 0
do nga=1,nuga
  do ngb=1,nga
    do ngc=1,ngb
      nla = nla+1
    end do
  end do
end do

nlb = 0
do nga=1,nuga
  do ngb=1,nga
    do ngc=1,nuga
      nlb = nlb+1
    end do
  end do
end do

return

end subroutine check_loops
