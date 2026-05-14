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

subroutine negp(ipdia,ipSigma,rout)

use ipPage, only: ipin, ipout, opout, W
use MCLR_Data, only: LuCIV, SS
use input_mclr, only: lRoots, nConf
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ipdia, ipSigma
real(kind=wp), intent(inout) :: rout(nConf)
integer(kind=iwp) :: i, iDisk
real(kind=wp), allocatable :: Tmp(:), Tmp2(:,:), Tmp3(:,:)
real(kind=wp), external :: DDot_

idisk = 0
call opout(ipdia)

call mma_allocate(Tmp,nConf,Label='Tmp')
call mma_allocate(Tmp2,2,lRoots,Label='Tmp2')
call mma_allocate(Tmp3,2,lRoots,Label='Tmp3')

call ipin(ipSigma)
do i=1,lroots
  call dDAFILE(luciv,2,Tmp,nconf,idisk)
  Tmp2(1,i) = DDOT_(nconf,rout,1,Tmp,1)
  Tmp2(2,i) = DDOT_(nconf,W(ipSigma)%A,1,Tmp,1)
end do
call ipout(ipsigma)
call dGeMV_('N',2*lroots,2*lroots,One,SS,2*lroots,Tmp2,1,Zero,Tmp3,1)

idisk = 0
call ipin(ipdia)
do i=1,lroots
  call dDAFILE(luciv,2,Tmp,nconf,idisk)
  call Exphinvv(W(ipdia)%A,Tmp,rout,One,Tmp3(1,i))
  rout(:) = rout(:)+Tmp3(2,i)*Tmp(:)
end do

call mma_deallocate(Tmp)
call mma_deallocate(Tmp2)
call mma_deallocate(Tmp3)

end subroutine negp
