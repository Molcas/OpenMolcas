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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine PutCMSFockOcc(FOccMO,nTri)

use MCLR_Data, only: ipMat, nDens
use input_mclr, only: nBas, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: FOccMO(nDens)
integer(kind=iwp), intent(in) :: nTri
integer(kind=iwp) :: iB, ijb, iS, jB
real(kind=wp), allocatable :: F(:), T(:)

call mma_allocate(F,nDens)
call mma_allocate(T,nDens)

F(:) = Zero
call Get_dArray_chk('FockOcc',F,nTri)
! WF Part
T(:) = FOccMO(:)
call TCMO(T,1,-2)
ijb = 0
do iS=1,nSym
  do ib=1,nbas(is)
    do jb=1,ib-1
      ijb = ijb+1
      F(ijb) = F(ijb)+T(ipmat(is,is)+nbas(is)*(jb-1)+ib-1)+T(ipmat(is,is)+nbas(is)*(ib-1)+jb-1)
    end do
    ijb = ijb+1
    F(ijb) = F(ijb)+T(ipmat(is,is)+nbas(is)*(ib-1)+ib-1)
  end do
end do
call Put_dArray('FockOcc',F,nDens)
call mma_deallocate(F)
call mma_deallocate(T)

end subroutine PutCMSFockOcc
