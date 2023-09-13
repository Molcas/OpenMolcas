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

logical function isitanint_cvb(a)

implicit real*8(a-h,o-z)
parameter(nallowed=12)
character*(*) a
character*(1) allowedchars(nallowed)
data allowedchars/'+','-','0','1','2','3','4','5','6','7','8','9'/

do ich=1,len_trim_cvb(a)
  do j=1,nallowed
    if (a(ich:ich) == allowedchars(j)) goto 100
  end do
  isitanint_cvb = .false.
  return
100 continue
end do
isitanint_cvb = .true.

return

end function isitanint_cvb
