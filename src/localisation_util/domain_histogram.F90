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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Domain_Histogram(iDomain,nAtom,nOcc,Title)
! Thomas Bondo Pedersen, January 2006.
!
! Purpose: print histogram of domain sizes.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtom, nOcc, iDomain(0:nAtom,nOcc)
character(len=*), intent(in) :: Title
integer(kind=iwp) :: i, i_max, i_min, iC, nC
real(kind=wp) :: Fac, Pct, x_ave
integer(kind=iwp), allocatable :: iCount(:)

if ((nAtom < 1) .or. (nOcc < 1)) return

i_min = iDomain(0,1)
i_max = iDomain(0,1)
x_ave = real(iDomain(0,1),kind=wp)
do i=2,nOcc
  i_min = min(i_min,iDomain(0,i))
  i_max = max(i_max,iDomain(0,i))
  x_ave = x_ave+real(iDomain(0,i),kind=wp)
end do
x_ave = x_ave/real(nOcc,kind=wp)

nC = i_max-i_min+1
call mma_allocate(iCount,nC,label='Dm_Histo')
iCount(:) = 0

call Cho_Head(Title,'=',80,u6)
write(u6,'(/,A,3X,I10,/,A,3X,I10,/,A,F13.2)') 'Minimum size:',i_min,'Maximum size:',i_max,'Average size:',x_ave

do i=1,nOcc
  iC = iDomain(0,i)-i_min+1
  iCount(iC) = iCount(iC)+1
end do

Fac = 1.0e2_wp/real(nOcc,kind=wp)
write(u6,*)
do iC=1,nC
  Pct = Fac*real(iCount(iC),kind=wp)
  write(u6,'(A,I10,A,I10,3X,F7.2,A)') 'Number with size',i_min+iC,':',iCount(iC),Pct,'%'
end do

call mma_deallocate(iCount)

end subroutine Domain_Histogram
