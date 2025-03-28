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

use ipPage, only: W
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two
use MCLR_Data, only: nConf1, ipCI
use MCLR_procedures, only: CISigma_sa
use input_mclr, only: nRoots, nCSF, State_Sym

implicit none
integer nh1, nh2, ipS1
real*8 h1(nh1), h2(nh2)
real*8, allocatable :: R(:,:)
real*8 rDum(1)
integer i, j
real*8, external :: DDot_

call CISigma_sa(0,state_sym,state_sym,h1,nh1,h2,nh2,rdum,1,ipCI,ipS1,.true.)

call ipin(ipS1)
call ipin(ipCI)

call DSCAL_(nroots*ncsf(STATE_SYM),Two,W(ipS1)%Vec,1)
call mma_allocate(R,[0,nroots-1],[0,nroots-1],label='R')

do i=0,nroots-1
  do j=0,nroots-1
    R(i,j) = ddot_(nconf1,W(ipS1)%Vec(1+nconf1*i),1,W(ipCI)%Vec(1+nconf1*j),1)
  end do
end do

do i=0,nroots-1
  do j=0,nroots-1
    call daxpy_(nconf1,-R(i,j),W(ipCI)%Vec(1+i*nconf1),1,W(ipS1)%Vec(1+j*nconf1),1)
  end do
end do

call mma_deallocate(R)

end subroutine Kap_CI
