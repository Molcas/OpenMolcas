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

subroutine mkp1(nEX,lst,rMat,rdiag)

use Constants, only: Zero
use negpre, only: LuCIV, P1
use stdalloc, only: mma_allocate, mma_deallocate
use input_mclr, only: lRoots, nConf, ERASSCF

implicit none
integer nEX
integer lst(nex)
real*8 rMat(*), rdiag(*)
real*8, allocatable :: Tmp1(:), Tmp2(:)
integer i, j, k, l, itri, idisk, jDisk, kk, ll
real*8 rtmp
! Statement functions
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

call mma_allocate(TMP1,nconf,Label='Tmp1')
call mma_allocate(TMP2,nconf,Label='Tmp2')

idisk = 0
do i=1,lroots
  call dDaFile(LuCIV,2,Tmp1,nconf,iDisk)
  jdisk = 0
  do j=1,i
    call dDafile(luciv,2,Tmp2,nconf,jDisk)
    rTmp = Zero
    do k=1,nex
      do l=1,nex
        kk = lst(k)
        ll = lst(l)
        rtmp = rtmp+Tmp1(kk)*Tmp2(ll)*rmat(itri(k,l))
      end do
    end do
    do k=1,nconf
      rtmp = rtmp+Tmp1(k)*Tmp2(k)*rdiag(k)
    end do
    if (i == j) rtmp = rtmp-ERASSCF(1)
    do k=1,nEx
      kk = lst(k)
      rtmp = rtmp-Tmp1(kk)*Tmp2(kk)*(rdiag(kk+1)-ERASSCF(1))
    end do
    P1(itri(i,j)) = rtmp
  end do
end do

call mma_deallocate(TMP2)
call mma_deallocate(TMP1)

end subroutine mkp1
