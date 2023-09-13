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

subroutine findmn_cvb(vec,n,vmn,imn)

implicit real*8(a-h,o-z)
dimension vec(n)

if (n > 0) then
  imn = 1
  vmn = vec(1)
  do i=2,n
    if (vec(i) < vmn) then
      imn = i
      vmn = vec(i)
    end if
  end do
else
  imn = 0
  vmn = 1d20
end if

return

end subroutine findmn_cvb
