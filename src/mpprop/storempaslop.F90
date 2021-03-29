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

subroutine StoreMpAsLop(nAtoms,ip_ANr,nB,ipT,ipTi,ipMP,lMax,ip_EC)

use MPProp_globals, only: Cor, iAtBoMltPlAd
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, ip_ANr, nB, ipT, ipTi, ipMP, lMax, ip_EC
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iAt1, iAt2, iAtK, iMu, ix, iy, j, kaunter, kompost, l, nSize1

!-- Let's fix the ip_ANr.

call Get_iArray('LP_A',iWork(ip_ANr),nAtoms)

!-- Let's fix the uber-simple T and T(-1).

kaunter = 0
do i=1,nB
  do j=1,nB
    Work(ipT+kaunter) = Zero
    Work(ipTi+kaunter) = Zero
    if (i == j) Work(ipT+kaunter) = One
    if (i == j) Work(ipTi+kaunter) = One
    kaunter = kaunter+1
  end do
end do

!-- Let's fix the expansion centres.

kaunter = 0
do i=1,nAtoms
  do j=1,i
    Work(ip_EC+kaunter*3+0) = Cor(1,i,j)
    Work(ip_EC+kaunter*3+1) = Cor(2,i,j)
    Work(ip_EC+kaunter*3+2) = Cor(3,i,j)
    kaunter = kaunter+1
  end do
end do

! Let's fix the multipole moments. Unlike LoProp, MpProp has here
! included the nuclei contribution, which we have to remove pronto
! to be compatible.

nSize1 = nAtoms*(nAtoms+1)/2
iMu = -1
do l=0,lMax
  kompost = 0
  do ix=l,0,-1
    do iy=l-ix,0,-1
      kompost = kompost+1
      iMu = iMu+1
      iAtK = 0
      do iAt1=1,nAtoms
        do iAt2=1,iAt1
          Work(ipMP+iAtK+nSize1*iMu) = Work(iAtBoMltPlAd(l)+nSize1*(kompost-1)+iAtK)
          iAtK = iAtK+1
        end do
        if (l == 0) then
          Work(ipMP+iAtK-1+nSize1*iMu) = Work(ipMP+iAtK-1+nSize1*iMu)-real(iWork(ip_ANr+iAt1-1),kind=wp)
        end if
      end do
    end do
  end do
end do

return

end subroutine StoreMpAsLop
