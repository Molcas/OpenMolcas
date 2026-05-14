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

subroutine mkguess_cvb()

use casvb_global, only: cvb, cvbdet, endvar, iapr, ipr, ixapr, kbasis, kbasiscvb, mxaobf, nbas_mo, norb, nvb, nvbinp, orbs, &
                        recinp, strtvb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: idum(1), ierr, ii, ioffs, iorb, iorb_ao, norb_ao
real(kind=wp) :: c, cnrm, dum(1)
integer(kind=iwp), allocatable :: irdorbs(:), itmp(:)
real(kind=wp), allocatable :: orbsao(:,:), tmp(:), tmp2(:,:)
real(kind=wp), parameter :: thresh = 1.0e-10_wp
real(kind=wp), external :: detm_cvb, dnrm2_, rand_cvb
logical(kind=iwp), external :: tstfile_cvb, & ! ... Files/Hamiltonian available ...
                               up2date_cvb    ! ... Make: up to date? ...

call mma_allocate(orbsao,nbas_mo,norb,label='orbsao')
call mma_allocate(irdorbs,norb,label='irdorbs')
irdorbs(:) = 0
! -- transfer from orbs if applicable -
! (Newly assigned memory => orbs will be zero)
do iorb=1,norb
  if (dnrm2_(norb,orbs(:,iorb),1) > thresh) then
    irdorbs(iorb) = 1
    orbsao(1:norb,iorb) = orbs(:,iorb)
  end if
end do
! -- restore from previous optim --
if (.not. up2date_cvb('RESTGS')) then
  if (up2date_cvb('WRITEGS')) call mkrestgs_cvb(orbsao,irdorbs,cvb,cvbdet,iapr,ixapr)
  call untouch_cvb('RESTGS')
end if
! -- read from file --
if (.not. up2date_cvb('STRTGS')) then
  call setstrtvb_cvb(strtvb)
  if (tstfile_cvb(strtvb)) call mkstrtgs_cvb(orbsao,irdorbs,cvb,strtvb)
  call untouch_cvb('STRTGS')
end if
! -- input --
if (.not. up2date_cvb('INPGS')) then
  call mma_allocate(itmp,norb,label='itmp')
  call rdioff_cvb(6,recinp,ioffs)
  call rdi_cvb(itmp,norb,recinp,ioffs)
  call rdioff_cvb(5,recinp,ioffs)
  do iorb=1,norb
    if (itmp(iorb) == 1) then
      ! MO basis ...
      irdorbs(iorb) = 1
      call rdlow_cvb(orbsao(1,iorb),norb,recinp,ioffs)
    else if (itmp(iorb) == 2) then
      ! AO basis ...
      irdorbs(iorb) = 2
      call rdlow_cvb(orbsao(1,iorb),nbas_mo,recinp,ioffs)
    end if
    ioffs = ioffs+mxaobf
  end do
  call mma_deallocate(itmp)

  call mma_allocate(tmp,nvbinp,label='tmp')
  call rdioff_cvb(7,recinp,ioffs)
  call rdrs_cvb(tmp,nvbinp,recinp,ioffs)
  if (dnrm2_(nvbinp,tmp,1) > thresh) then
    call rdioff_cvb(3,recinp,ioffs)
    call rdis_cvb(idum,1,recinp,ioffs)
    kbasiscvb = idum(1)
    cvb(1:nvbinp) = tmp(:)
  end if
  call mma_deallocate(tmp)

  call untouch_cvb('INPGS')
end if
! -- semi-random --
! Leading diagonal, random but positive orbital overlaps:
dum = rand_cvb(0.777_wp)
c = 0.1_wp
do iorb=1,norb
  if (irdorbs(iorb) == 0) then
    irdorbs(iorb) = 1
    do ii=1,norb
      orbsao(ii,iorb) = c*rand_cvb(Zero)
      if (ii == iorb) orbsao(ii,iorb) = One
    end do
  else
    ! Dummy calls to RAND to get consistent guesses:
    do ii=1,norb
      dum = rand_cvb(Zero)
    end do
  end if
end do

! Collect orbitals and transform AO -> MO:
norb_ao = 0
do iorb=1,norb
  if (irdorbs(iorb) == 1) then
    orbs(:,iorb) = orbsao(1:norb,iorb)
  else if (irdorbs(iorb) == 2) then
    norb_ao = norb_ao+1
    if (norb_ao /= iorb) orbsao(:,norb_ao) = orbsao(:,iorb)
  end if
end do
call mma_allocate(tmp2,norb,norb_ao,label='tmp2')
call ao2mo_cvb(orbsao,tmp2,norb_ao)
iorb_ao = 0
do iorb=1,norb
  if (irdorbs(iorb) == 2) then
    iorb_ao = iorb_ao+1
    orbs(:,iorb) = tmp2(:,iorb_ao)
  end if
end do
call mma_deallocate(tmp2)
call mma_deallocate(orbsao)
call mma_deallocate(irdorbs)

ierr = 0
call nize_cvb(orbs,norb,dum,norb,0,ierr)

if (abs(detm_cvb(orbs,norb)) < 1.0e-8_wp) then
  dum = rand_cvb(0.777_wp)
  c = 0.1_wp
  do iorb=1,norb
    do ii=1,norb
      orbs(ii,iorb) = orbs(ii,iorb)+c*(One-Two*rand_cvb(Zero))
    end do
  end do
  if (abs(detm_cvb(orbs,norb)) < 1.0e-8_wp) then
    if (ipr(1) >= 0) write(u6,'(a)') ' Starting orbital guess was near-singular - using semi-random guess instead.'
    dum = rand_cvb(0.777_wp)
    c = 0.1_wp
    do iorb=1,norb
      do ii=1,norb
        orbs(ii,iorb) = c*rand_cvb(Zero)
        if (ii == iorb) orbs(ii,iorb) = One
      end do
    end do
  else
    if (ipr(1) >= 0) write(u6,'(a)') ' Starting orbital guess was near-singular - scrambling orbital coefficients.'
  end if
  ierr = 0
  call nize_cvb(orbs,norb,dum,norb,0,ierr)
end if

ierr = 0
call nize_cvb(orbs,norb,dum,norb,0,ierr)

! Perfect-pairing spin function(s):
if (dnrm2_(nvb,cvb,1) < thresh) then
  kbasiscvb = kbasis
  call ppgs_cvb(cvb)
end if

cnrm = dnrm2_(nvb,cvb,1)
if (cnrm < thresh) then
  write(u6,*) ' Fatal error - starting structure coefficients all zero !'
  call abend_cvb()
end if

if (kbasiscvb /= kbasis) then
  call mktrnspn_cvb()
  call untouch_cvb('TRNSPN')
end if

!if (ploc) call rtransf_plc(orbs,cvb)

if ((ipr(1) >= 2) .and. (.not. endvar)) then
  write(u6,'(/,a)') ' Wavefunction guess :'
  call report_cvb(orbs,norb)
  write(u6,'(/,a)') ' Structure coefficients :'
  write(u6,'(a)') ' ------------------------'
  call vecprint_cvb(cvb,nvb)
end if

return

end subroutine mkguess_cvb
