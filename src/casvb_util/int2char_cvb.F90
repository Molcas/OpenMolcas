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

subroutine int2char_cvb(a,intx,iform)
! Simulates internal write:    write(a,'(i<iform>)') intx

use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(out) :: a
integer(kind=iwp), intent(in) :: intx, iform
real(kind=wp) :: dum
integer(kind=iwp) :: i, ia, iamax, idum, int2, la, numb
character, parameter :: blnk = ' ', cnumb(0:9) = ['0','1','2','3','4','5','6','7','8','9'], minus = '-'

la = len(a)
if (iform > la) then
  write(u6,*) ' Illegal call to int2char_cvb:',iform,la
  call abend_cvb()
end if
dum = log10(real(max(abs(intx),1),kind=wp))
idum = nint(dum)
if (abs(intx) >= 10**idum) then
  iamax = idum+1
else
  iamax = idum
end if
if (intx < 0) iamax = iamax+1
if (iamax > iform) then
  write(u6,*) ' Integer too large in int2char_cvb:',intx,iform
  call abend_cvb()
end if
a(1:iform-iamax) = ''
ia = iform-iamax
if (intx < 0) then
  ia = ia+1
  a(ia:ia) = minus
  iamax = iamax-1
end if
int2 = abs(intx)
do i=iamax-1,0,-1
  ia = ia+1
  numb = int2/(10**i)
  a(ia:ia) = cnumb(numb)
  int2 = int2-numb*(10**i)
end do
if (intx == 0) a(iform:iform) = cnumb(0)

return

end subroutine int2char_cvb
