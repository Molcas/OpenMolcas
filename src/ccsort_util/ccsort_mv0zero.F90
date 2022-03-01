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

subroutine ccsort_mv0zero(DD,LENGTH,MAT)

integer DD
integer LENGTH
real*8 MAT(1:DD)
integer INIT
real*8 ZERO
data ZERO/0.0D+00/

! ...loop over all elements

do INIT=1,LENGTH
  MAT(INIT) = ZERO
end do

return

end subroutine ccsort_mv0zero
