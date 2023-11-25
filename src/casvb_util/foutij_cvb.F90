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

subroutine foutij_cvb(fij,ni,nj,a1,a2)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ni, nj
real(kind=wp), intent(in) :: fij(ni,nj)
character(len=*), intent(in) :: a1, a2
character(len=46) :: b2
character(len=15) :: b1
character(len=12) :: b3
integer(kind=iwp) :: i, ichar0, ichar9, ipos, j, jpos
logical(kind=iwp) :: done
real(kind=wp), parameter :: hge = 1.0e20_wp

b1 = a1
b2 = a2
b3 = '     ......'
write(u6,'(/,1x,3a)') b1,b2,b3
b2 = ' '
! Find IPOS/JPOS : position of I/J indices in string
ichar0 = ichar('0')
ichar9 = ichar('9')
done = .false.
do jpos=15,1,-1
  if ((ichar(b1(jpos:jpos)) >= ichar0) .and. (ichar(b1(jpos:jpos)) <= ichar9)) then
    done = .true.
    exit
  end if
end do
if (.not. done) then
  write(u6,*) ' Fatal error in FOUTIJ!'
  call abend_cvb()
end if
done = .false.
do ipos=jpos-1,1,-1
  if ((ichar(b1(ipos:ipos)) >= ichar0) .and. (ichar(b1(ipos:ipos)) <= ichar9)) then
    done = .true.
    exit
  end if
end do
if (.not. done) then
  write(u6,*) ' Fatal error in FOUTIJ!'
  call abend_cvb()
end if
do j=1,nj
  do i=1,ni
    if (abs(fij(i,j)) /= hge) then
      write(b1(ipos:ipos),'(i1)') i
      write(b1(jpos:jpos),'(i1)') j
      write(b3,'(es12.4)') fij(i,j)
      write(u6,'(1x,3a)') b1,b2,b3
    end if
  end do
end do

return

end subroutine foutij_cvb
