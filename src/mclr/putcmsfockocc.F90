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

use MCLR_Data, only: nDens2, ipMat
use input_mclr, only: nSym, nBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One

implicit none
! Input:
integer nTri
real*8, dimension(nDens2) :: FOccMO
! Auxiliaries
real*8, dimension(:), allocatable :: F, T, F_n
integer ijb, iS, iB, jB

call mma_allocate(F,nDens2)
call mma_allocate(T,nDens2)
call mma_allocate(F_n,nDens2)

call FZero(F,nDens2)
call FZero(F_n,nDens2)
call Get_dArray_chk('FockOcc',F,nTri)
! WF Part
call DCopy_(nDens2,FOccMO,1,T,1)
call TCMO(T,1,-2)
ijb = 0
do iS=1,nSym
  do ib=1,nbas(is)
    do jb=1,ib-1
      ijb = ijb+1
      F_n(ijb) = T(ipmat(is,is)+nbas(is)*(JB-1)+IB-1)+T(ipmat(is,is)+nbas(is)*(IB-1)+JB-1)
    end do
    ijb = ijb+1
    F_n(ijb) = T(ipmat(is,is)+nbas(is)*(iB-1)+IB-1)
  end do
end do
call Daxpy_(nDens2,One,F_n,1,F,1)
call Put_dArray('FockOcc',F,nDens2)
call mma_deallocate(F)
call mma_deallocate(T)
call mma_deallocate(F_n)

end subroutine PutCMSFockOcc
