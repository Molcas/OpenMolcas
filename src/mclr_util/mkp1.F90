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

use Index_Functions, only: iTri
use MCLR_Data, only: LuCIV, P1
use input_mclr, only: ERASSCF, lRoots, nConf
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nEX, lst(nex)
real(kind=wp), intent(in) :: rMat(*), rdiag(*)
integer(kind=iwp) :: i, idisk, j, jDisk, k, kk, l
real(kind=wp) :: rtmp
real(kind=wp), allocatable :: Tmp1(:), Tmp2(:)

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
        rtmp = rtmp+Tmp1(lst(k))*Tmp2(lst(l))*rmat(iTri(k,l))
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
    P1(iTri(i,j)) = rtmp
  end do
end do

call mma_deallocate(TMP2)
call mma_deallocate(TMP1)

end subroutine mkp1
