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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine icomb_cvb(i1,i2,icomb_cvbval)

implicit real*8(a-h,o-z)
save one, half
data one/1.d0/,half/0.5d0/

icomb_cvbval = 0

! Special cases - return ICOMB_CVB=0:
if ((i1 < 0) .or. (i2 < 0) .or. (i1 < i2)) then
  icomb_cvbval = 0
  return
end if
! I3 is I2 but always less than I1/2:
i3 = (i1-abs(i1-2*i2))/2
comb = one
do j=1,i3
  comb = comb/dble(j)
  comb = comb*dble(i1-j+1)
end do
icomb_cvbval = nint(comb)
! If integer overflow - return ICOMB_CVB=-1:
if (abs(dble(icomb_cvbval)-comb) > half) icomb_cvbval = -1

return

end subroutine icomb_cvb
