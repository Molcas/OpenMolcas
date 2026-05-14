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

subroutine getmo2_cvb(cmo,cmo2,ic)

use casvb_global, only: iact_mo, nact_mo, nbas_mo, nbasf_mo, nbasi_mo, nbasisq_mo, nbassqf_mo, nsym_mo
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: cmo(nbas_mo,nbas_mo)
real(kind=wp), intent(inout) :: cmo2(nbas_mo,nbas_mo)
integer(kind=iwp), intent(in) :: ic
integer(kind=iwp) :: iorb, isk, jbas
real(kind=wp), allocatable :: cmoblk(:)

! Construct full matrix of MOs in symmetry-adapted AO basis:
call mma_allocate(cmoblk,nbasisq_mo,label='cmoblk')
call getmoblk_cvb(cmoblk)
cmo(:,:) = Zero
do isk=1,nsym_mo
  do jbas=1,nbasi_mo(isk)
    cmo(nbasf_mo(isk)+1:nbasf_mo(isk)+nbasi_mo(isk),jbas+nbasf_mo(isk)) = &
      cmoblk(nbassqf_mo(isk)+(jbas-1)*nbasi_mo(isk)+1:nbassqf_mo(isk)+jbas*nbasi_mo(isk))
  end do
end do
call mma_deallocate(cmoblk)

if (mod(ic,2) == 1) then
  call mxinv_cvb(cmo,nbas_mo)
  call dgetmi(cmo,nbas_mo,nbas_mo)
end if

if (ic >= 2) then
  do iorb=1,nact_mo
    cmo2(:,iorb) = cmo(:,iact_mo(iorb))
  end do
end if

return

end subroutine getmo2_cvb
