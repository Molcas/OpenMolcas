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

subroutine gentkin(L,TKIN,nprims,exponents,rootOVLPinv)
!bs subroutine to generate the kinetic energy

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: L, nprims
#include "para.fh"
real(kind=wp) :: TKIN(nprims,nprims), exponents(*), rootOVLPinv(MxprimL,MxprimL)
integer(kind=iwp) :: irun, irun1, irun2, jrun, krun
real(kind=wp), allocatable :: dummy(:,:), dummy2(:,:)
real(kind=wp), external :: Tkinet

call mma_allocate(dummy,MxprimL,MxprimL,label='dummy')
call mma_allocate(dummy2,MxprimL,MxprimL,label='dummy2')

!bs one triangular part of the matrix
do irun2=1,nprims
  do irun1=1,irun2
    dummy(irun1,irun2) = Tkinet(l,exponents(irun1),exponents(irun2))
  end do
end do
!bs copy to the other triangular part....
do irun2=1,nprims-1
  do irun1=irun2+1,nprims
    dummy(irun1,irun2) = dummy(irun2,irun1)
  end do
end do
!bs now transform by rootovlp*dummy*rootovlp
TKIN(1:nprims,1:nprims) = Zero
dummy2(1:nprims,1:nprims) = Zero
do irun=1,nprims
  do jrun=1,nprims
    do krun=1,nprims
      dummy2(irun,jrun) = dummy2(irun,jrun)+dummy(irun,krun)*rootovlpinv(krun,jrun)
    end do
  end do
end do
do irun=1,nprims
  do jrun=1,nprims
    do krun=1,nprims
      Tkin(irun,jrun) = Tkin(irun,jrun)+dummy2(krun,jrun)*rootovlpinv(irun,krun)
    end do
  end do
end do

call mma_deallocate(dummy)
call mma_deallocate(dummy2)

return

end subroutine gentkin
