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

subroutine mkp1inv(rdia)

use Index_Functions, only: iTri
use MCLR_Data, only: LuCIV, P1Inv
use input_mclr, only: lRoots, nConf
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: rdia(*)
integer(kind=iwp) :: i, iDisk, j, jDisk
real(kind=wp), allocatable :: TMP1(:), TMP2(:)
real(kind=wp), external :: DDot_

call mma_allocate(TMP1,nconf,Label='TMP1')
call mma_allocate(TMP2,nconf,Label='TMP2')

idisk = 0
do i=1,lroots
  jdisk = idisk
  call dDaFile(LuCIV,2,Tmp1,nconf,iDisk)
  call ExpHinvv(rdia,Tmp1,Tmp1,Zero,One)
  do j=i,lroots
    call dDafile(luciv,2,Tmp2,nconf,jDisk)
    p1INV(iTri(i,j)) = DDOT_(nconf,Tmp2,1,Tmp1,1)
  end do
end do

call mma_deallocate(TMP1)
call mma_deallocate(TMP2)

end subroutine mkp1inv
