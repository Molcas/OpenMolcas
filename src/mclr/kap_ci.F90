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

subroutine Kap_CI(h1,nh1,h2,nh2,ipS1)

use ipPage, only: ipin, W
use MCLR_Data, only: ipCI, nConf1
use MCLR_procedures, only: CISigma_sa
use input_mclr, only: nCSF, nRoots, State_Sym
use ISRotation, only: unequal_SA
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nh1, nh2, ipS1
real(kind=wp), intent(in) :: h1(nh1), h2(nh2)
integer(kind=iwp) :: i, j
real(kind=wp) :: rDum(1)
real(kind=wp), allocatable :: R(:,:)
real(kind=wp), external :: DDot_

call CISigma_sa(0,state_sym,state_sym,h1,nh1,h2,nh2,rdum,1,ipCI,ipS1,.true.)

call ipin(ipS1)
call ipin(ipCI)

W(ipS1)%A(1:nroots*ncsf(STATE_SYM)) = Two*W(ipS1)%A(1:nroots*ncsf(STATE_SYM))

! if not equally state-averaged, do not project out the internal space contributions here
if (unequal_SA) return

call mma_allocate(R,nroots,nroots,label='R')

do i=1,nroots
  do j=1,nroots
    R(i,j) = ddot_(nconf1,W(ipS1)%A(1+nconf1*(i-1)),1,W(ipCI)%A(1+nconf1*(j-1)),1)
  end do
end do

do i=1,nroots
  do j=1,nroots
    W(ipS1)%A((j-1)*nconf1+1:j*nconf1) = W(ipS1)%A((j-1)*nconf1+1:j*nconf1)-R(i,j)*W(ipCI)%A((i-1)*nconf1+1:i*nconf1)
  end do
end do

call mma_deallocate(R)

end subroutine Kap_CI
