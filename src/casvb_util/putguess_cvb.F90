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

subroutine putguess_cvb(orbs,cvb,recn)

use casvb_global, only: endvar, ifmos, ipr, kbasiscvb, nbas_mo, norb, nvb, ploc, variat
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

#include "intent.fh"
implicit none
real(kind=wp), intent(_IN_) :: orbs(norb,*), cvb(*)
real(kind=wp), intent(in) :: recn
integer(kind=iwp) :: i, ierr, ioffs_cvb, ioffs_orbs, ioffs_orbsao, ioffs_orbslao, iorb, kbasiscvb1, nbas_mo1, norb1, nvb1
logical(kind=iwp) :: use_ao
real(kind=wp), allocatable :: a(:,:), b(:,:), c(:), orbsao(:,:)
real(kind=wp), external :: dnrm2_

call wrheader_cvb(recn,norb,nbas_mo,nvb,kbasiscvb,ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
call rdheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
do iorb=1,norb
  call wrgspr_cvb(recn,orbs(:,iorb),iorb,norb,1,ierr)
end do
call wrgspr_cvb(recn,cvb,1,nvb,2,ierr)
use_ao = ifmos .and. ((.not. variat) .or. (variat .and. endvar))
if (use_ao) then
  call mma_allocate(orbsao,nbas_mo,norb)
  call mo2ao_cvb(orbs,orbsao,norb)
  do iorb=1,norb
    call wrgspr_cvb(recn,orbsao(:,iorb),iorb,nbas_mo,3,ierr)
  end do
  if (ipr(5) >= 2) then
    write(u6,'(/,a)') ' VB orbitals in AO basis :'
    write(u6,'(a)') ' -------------------------'
    call mxprint_cvb(orbsao,nbas_mo,norb,0)
  end if
  if (ploc) then
    call untested('putguess_cvb: ploc')
    call mma_allocate(a,norb,norb,label='a')
    call mma_allocate(b,norb,norb,label='b')
    call mma_allocate(c,norb,label='c')
    !call getr_plc(a)
    call dgetmi(a,norb,norb)
    call mxatb_cvb(a,orbs,norb,norb,norb,b)
    call lmo2ao_cvb(b,orbsao,norb)
    do iorb=1,norb
      call wrgspr_cvb(recn,orbsao(:,iorb),iorb,nbas_mo,4,ierr)
    end do
    if (ipr(5) >= 2) then
      write(u6,'(/,a)') ' Original localized VB orbitals in AO basis :'
      write(u6,'(a)') ' --------------------------------------------'
      call mxprint_cvb(orbsao,nbas_mo,norb,0)
    end if
    do i=1,norb
      c(i) = dnrm2_(norb,b(:,i),1)
      b(:,i) = b(:,i)/c(i)
    end do
    if (ipr(5) >= 2) then
      write(u6,'(/,a)') ' Norms of original localized VB orbitals :'
      write(u6,'(a)') ' -----------------------------------------'
      call mxprint_cvb(c,1,norb,0)
    end if
    call mma_deallocate(a)
    call mma_deallocate(b)
    call mma_deallocate(c)
  end if
  call mma_deallocate(orbsao)
end if

return

end subroutine putguess_cvb
