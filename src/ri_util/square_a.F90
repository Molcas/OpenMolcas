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

subroutine Square_A(Lu,nB,MaxMem_,Force_out_of_Core)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lu, nB, MaxMem_
logical(kind=iwp), intent(in) :: Force_out_of_Core
integer(kind=iwp) :: iAddr, iAddr1, iAddr2, iAddrs, iB, Inc, jB, kB, MaxMem, mB, nBuff, nMem
real(kind=wp), allocatable :: Buf(:,:)

if (nB == 0) return
MaxMem = MaxMem_
nMem = nB**2
if (Force_Out_of_Core) MaxMem = nMem/3

if (nMem <= MaxMem) then

  ! In-core case

  call mma_allocate(Buf,nMem,1,Label='Buf')
  iAddr = 0
  call dDaFile(Lu,2,Buf(:,1),nMem,iAddr)
  call In_place_Square(Buf(:,1),nB)
  iAddr = 0
  call dDaFile(Lu,1,Buf(:,1),nMem,iAddr)

else

  ! Out-of-core case

  nBuff = MaxMem/2
  call mma_allocate(Buf,nBuff,2,Label='Buf')

  Inc = nBuff/nB
  iAddr1 = 0
  do iB=1,nB,Inc
    mB = min(Inc,nB-iB+1)
    iAddrs = iAddr1
    call dDaFile(Lu,2,Buf(:,1),nB*mB,iAddr1)

    iAddr2 = iAddr1
    do jB=iB,nB,Inc
      kB = min(Inc,nB-jB+1)

      if (jB == iB) then
        call In_place_Diag(Buf(:,1),nB,iB,iB+mB-1)
      else
        call dDaFile(Lu,2,Buf(:,2),nB*kB,iAddr2)
        call Off_Diagonal(Buf(:,1),nB,iB,iB+mB-1,Buf(:,2),jB,jB+kB-1)
      end if
    end do

    iAddr1 = iAddrs
    call dDaFile(Lu,1,Buf(:,1),nB*mB,iAddr1)

  end do
end if
call mma_deallocate(Buf)

return

end subroutine Square_A
