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

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nAtom, nOcc, iDomain(0:nAtom,nOcc)
character(len=*) :: Title
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i_max, i_min, iC, ip_iCount, l_iCount
real(kind=wp) :: Fac, Pct, x_ave

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

l_iCount = i_max-i_min+1
call GetMem('Dm_Histo','Allo','Inte',ip_iCount,l_iCount)

call Cho_Head(Title,'=',80,u6)
write(u6,'(/,A,3X,I10,/,A,3X,I10,/,A,F13.2)') 'Minimum size:',i_min,'Maximum size:',i_max,'Average size:',x_ave
call Domain_Histo1(iDomain,nAtom,nOcc,iWork(ip_iCount),i_min,i_max)
Fac = 1.0e2_wp/real(nOcc,kind=wp)
Pct = Fac*real(iWork(ip_iCount),kind=wp)
i = i_min
write(u6,'(/,A,I10,A,I10,3X,F7.2,A)') 'Number with size',i,':',iWork(ip_iCount),Pct,'%'
do iC=2,l_iCount
  Pct = Fac*real(iWork(ip_iCount-1+iC),kind=wp)
  i = i+1
  write(u6,'(A,I10,A,I10,3X,F7.2,A)') 'Number with size',i,':',iWork(ip_iCount-1+iC),Pct,'%'
end do

call GetMem('Dm_Histo','Free','Inte',ip_iCount,l_iCount)

end subroutine Domain_Histogram
