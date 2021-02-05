!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

subroutine set_l_kriging(lv,nInter_In)

use kriging_mod, only: l, nInter
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nInter_In
real(kind=wp) :: lv(nInter_In)

! Set the characteristic length of all the components of the coordintes.

if (nInter_In == nInter) then
  l(:) = lv(:)
else if (nInter == 1) then
  l(:) = lv(1)
else
  write(u6,*) 'setlkriging: illegal nInter value.'
  call Abend()
end if

! Generate the covariance matrix

call covarMatrix()

! Form the inverse of the covariance matrix times the generalized value vector.

call kriging_model()

end subroutine set_l_kriging
