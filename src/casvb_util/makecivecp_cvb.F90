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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: civec(ndet), civecp(ndet), orbs(norb,norb)
real(kind=wp), allocatable :: gjorb(:), owrk(:,:)
integer(kind=iwp), external :: ihlf_cvb
logical(kind=iwp), external :: tstcnt_cvb ! ... Content of CI vectors ...

if (tstcnt_cvb(civecp,3)) return

call mma_allocate(owrk,norb,norb,label='owrk')
call mma_allocate(gjorb,norb*norb+ihlf_cvb(norb+2*norb*norb),label='gjorb')
call transp_cvb(orbs,owrk,norb,norb)
call gaussj_cvb(owrk,gjorb)
if (memplenty) then
  call getci_cvb(civec)
  call cicopy_cvb(civec,civecp)
else
  call cird_cvb(civecp,61001.2_wp)
end if
call applyt_cvb(civecp,gjorb)
call mma_deallocate(owrk)
call mma_deallocate(gjorb)

call setcnt_cvb(civecp,3)

return

end subroutine makecivecp_cvb
