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

subroutine fouti_cvb(fi,ni,a1,a2)

implicit real*8(a-h,o-z)
!logical l
character*(*) a1, a2
character*15 b1
character*46 b2
character*12 b3
dimension fi(ni)
save huge
data huge/1d20/

b1 = a1
b2 = a2
b3 = '     ......'
write(6,'(/,1x,3a)') b1,b2,b3
b2 = ' '
! Find IPOS : position of I index in string
ichar0 = ichar('0')
ichar9 = ichar('9')
do ipos=15,1,-1
  if ((ichar(b1(ipos:ipos)) >= ichar0) .and. (ichar(b1(ipos:ipos)) <= ichar9)) goto 100
end do
write(6,*) ' Fatal error in FOUTI!'
call abend_cvb()
100 continue
do i=1,ni
  if (abs(fi(i)) /= huge) then
    write(b1(ipos:ipos),'(i1)') i
    write(b3,'(d12.4)') fi(i)
    write(6,'(1x,3a)') b1,b2,b3
  end if
end do

return

end subroutine fouti_cvb
