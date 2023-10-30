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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine makecivecp_cvb(civec,civecp,orbs)
! Construct CIVECP:

use casvb_global, only: gjorb_type, icnt_ci, memplenty, ndet, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: civec(0:ndet), civecp(0:ndet)
real(kind=wp), intent(in) :: orbs(norb,norb)
integer(kind=iwp) :: icivecp
type(gjorb_type) :: gjorb
real(kind=wp), allocatable :: owrk(:,:)

icivecp = nint(civecp(0))
if (icnt_ci(icivecp) == 3) return

call mma_allocate(owrk,norb,norb,label='owrk')
call mma_allocate(gjorb%r,norb,norb,label='gjorb%r')
call mma_allocate(gjorb%i1,norb,label='gjorb%i1')
call mma_allocate(gjorb%i2,2,norb*norb,label='gjorb%i2')
call trnsps(norb,norb,orbs,owrk)
call gaussj_cvb(owrk,gjorb)
if (memplenty) then
  call getci_cvb(civec)
  call cicopy_cvb(civec,civecp)
else
  call cird_cvb(civecp,61001.2_wp)
end if
call applyt_cvb(civecp,gjorb)
call mma_deallocate(owrk)
call mma_deallocate(gjorb%r)
call mma_deallocate(gjorb%i1)
call mma_deallocate(gjorb%i2)

icnt_ci(icivecp) = 3

return

end subroutine makecivecp_cvb
