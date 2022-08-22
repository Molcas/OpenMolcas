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

subroutine WrHDsk(Hess,nGrad)

use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nGrad
real(kind=wp) :: Hess(nGrad*(nGrad+1)/2)
#include "Molcas.fh"
#include "disp.fh"
integer(kind=iwp) :: i, idum, iG, iG1, iG2, iGrad1, iGrad2, iIrrep, iOpt, ip_Acc, iRc, mH, nH
character(len=8) :: Label
real(kind=wp), allocatable :: EVal(:), EVec(:), HStat(:), Temp(:)

call mma_allocate(Temp,nGrad**2,Label='Temp')
nH = 0
do iIrrep=0,nIrrep-1
  nH = nH+lDisp(iIrrep)
end do
call mma_allocate(HStat,nH,Label='HStat')

! Reorder Hessian to lower triangular form

iGrad1 = 1
iGrad2 = 0
iG = 0
ip_Acc = 1
do iIrrep=0,nIrrep-1
  iGrad2 = iGrad2+lDisp(iIrrep)

  do iG1=iGrad1,iGrad2
    do iG2=iGrad1,iG1
      iG = iG+1
      Temp(iG) = Hess(iG1*(iG1-1)/2+IG2)
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Diagonalize and keep eigenvalues for check facility

  mH = lDisp(iIrrep)
  call mma_allocate(EVal,mH*(mH+1)/2,Label='EVal')
  call mma_allocate(EVec,mH*mH,Label='EVec')

  call dcopy_(mH*(mH+1)/2,Temp,1,EVal,1)
  call dcopy_(mH*mH,[Zero],0,EVec,1)
  call dcopy_(mH,[One],0,EVec,mH+1)

  ! Compute eigenvalues and eigenvectors

  call Jacob(EVal,EVec,mH,mH)
  call Jacord(EVal,EVec,mH,mH)

  do i=1,mH
    HStat(ip_Acc) = EVal(i*(i+1)/2)
    ip_Acc = ip_Acc+1
  end do

  call mma_deallocate(EVec)
  call mma_deallocate(EVal)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  iGrad1 = iGrad1+lDisp(iIrrep)
end do

! Write eigenvalues to the check file.

call Add_Info('HStat',HStat,nH,5)

iRc = -1
iOpt = 0
Label = 'StatHess'
call dWrMck(iRC,iOpt,Label,idum,Temp,idum)
if (iRc /= 0) then
  write(u6,*) 'WrHDsk: Error writing to MCKINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
call mma_deallocate(HStat)
call mma_deallocate(Temp)

return

end subroutine WrHDsk
