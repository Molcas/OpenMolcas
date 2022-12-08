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

subroutine CD_AInv(A,n,AInV,Thr_CD)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(_IN_) :: A(n,n)
real(kind=wp), intent(out) :: AInv(n,n)
real(kind=wp), intent(in) :: Thr_CD
integer(kind=iwp) :: iDisk, Lu_A, Lu_Q, m
real(kind=wp), allocatable :: ADiag(:), QVec(:,:)
#ifdef _ACCURACY_
real(kind=wp), allocatable :: Tmp(:,:), Tmp2(:,:)
#endif
integer(kind=iwp), external :: IsFreeUnit

call mma_allocate(ADiag,n,Label='ADiag')

Lu_A = IsFreeUnit(77)
call DaName_MF_WA(Lu_A,'AMat09')

iDisk = 0
call dDaFile(Lu_A,1,A,n**2,iDisk)

!call RecPrt('A',' ',A,n,n)

Lu_Q = IsFreeUnit(78)
call DaName_MF_WA(Lu_Q,'QMat09')

call dcopy_(n,A,n+1,ADiag,1)

call CD_AInv_Inner(n,m,ADiag,Lu_A,Lu_Q,Thr_CD)

call mma_deallocate(ADiag)

call mma_allocate(QVec,n,m,Label='QVec')

iDisk = 0
call dDaFile(Lu_Q,2,QVec,n*m,iDisk)

!call RecPrt('QVec','(6G20.10)',QVec,n,m)
call DGEMM_('N','T',n,n,m,One,QVec,n,QVec,n,Zero,AInv,n)
!call RecPrt('AInv',' ',AInv,n,n)
call DaEras(Lu_Q)
call mma_deallocate(QVec)
!                                                                      *
!***********************************************************************
!                                                                      *
! Check the accuracy I-AA^1

#ifdef _ACCURACY_
call mma_allocate(Tmp,n,n,Label='Tmp')
call mma_allocate(Tmp2,n,n,Label='Tmp2')

! I
call unitmat(Tmp,n)
Tmp2(:,:) = Tmp
! I-AA^-1
call DGEMM_('N','N',n,n,n,-One,A,n,AInv,n,One,Tmp,n)
call RecPrt('I-AA^-1','(6G20.12)',Tmp,n,n)

call DGEMM_('N','N',n,n,n,One,A,n,AInv,n,Zero,Tmp,n)

call DGEMM_('N','N',n,n,n,-One,Tmp,n,Tmp,n,One,Tmp2,n)
call RecPrt('I-AA^-1AA^-1','(6G20.12)',Tmp2,n,n)

call mma_deallocate(Tmp2)
call mma_deallocate(Tmp)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine CD_AInv
