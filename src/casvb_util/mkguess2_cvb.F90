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

subroutine mkguess2_cvb(orbs,cvb,irdorbs,orbsao)

implicit real*8(a-h,o-z)
! ... Files/Hamiltonian available ...
logical, external :: tstfile_cvb
! ... Make: up to date? ...
logical, external :: up2date_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
#include "mo_cvb.fh"
dimension orbs(norb,norb), cvb(nvb)
dimension irdorbs(norb), orbsao(nbas_mo,norb)
dimension idum(1), dum(1)
save thresh
data thresh/1d-10/

call izero(irdorbs,norb)
! -- transfer from orbs if applicable -
! (Newly assigned memory => orbs will be zero)
do iorb=1,norb
  if (dnrm2_(norb,orbs(1,iorb),1) > thresh) then
    irdorbs(iorb) = 1
    call fmove_cvb(orbs(1,iorb),orbsao(1,iorb),norb)
  end if
end do
! -- restore from previous optim --
if (.not. up2date_cvb('RESTGS')) then
  if (up2date_cvb('WRITEGS')) then
    call rdi_cvb(idum,1,recn_tmp04,0)
    ndetvb1 = idum(1)
    i1 = mstacki_cvb(ndetvb1)
    i2 = mstackr_cvb(ndetvb1)
    call mkrestgs_cvb(orbsao,irdorbs,cvb,work(lw(9)),iwork(ll(11)),iwork(ll(12)),iwork(i1),work(i2))
    call mfreei_cvb(i1)
  end if
  call untouch_cvb('RESTGS')
end if
! -- read from file --
if (.not. up2date_cvb('STRTGS')) then
  call setstrtvb_cvb(strtvb)
  if (tstfile_cvb(strtvb)) call mkstrtgs_cvb(orbsao,irdorbs,cvb,strtvb,kbasiscvb)
  call untouch_cvb('STRTGS')
end if
! -- input --
if (.not. up2date_cvb('INPGS')) then
  i1 = mstacki_cvb(norb)
  call rdioff_cvb(6,recinp,ioffs)
  call rdi_cvb(iwork(i1),norb,recinp,ioffs)
  call rdioff_cvb(5,recinp,ioffs)
  do iorb=1,norb
    if (iwork(iorb+i1-1) == 1) then
      ! MO basis ...
      irdorbs(iorb) = 1
      call rdr_cvb(orbsao(1,iorb),norb,recinp,ioffs)
    else if (iwork(iorb+i1-1) == 2) then
      ! AO basis ...
      irdorbs(iorb) = 2
      call rdr_cvb(orbsao(1,iorb),nbas_mo,recinp,ioffs)
    end if
    ioffs = ioffs+mxaobf
  end do
  call mfreei_cvb(i1)

  i1 = mstackr_cvb(nvbinp)
  call rdioff_cvb(7,recinp,ioffs)
  call rdrs_cvb(work(i1),nvbinp,recinp,ioffs)
  if (dnrm2_(nvbinp,work(i1),1) > thresh) then
    call rdioff_cvb(3,recinp,ioffs)
    call rdis_cvb(idum,1,recinp,ioffs)
    kbasiscvb = idum(1)
    call fmove_cvb(work(i1),work(lv(2)),nvbinp)
  end if
  call mfreer_cvb(i1)

  call untouch_cvb('INPGS')
end if
! -- semi-random --
! Leading diagonal, random but positive orbital overlaps:
dum = rand_cvb(.777d0)
c = 1d-1
do iorb=1,norb
  if (irdorbs(iorb) == 0) then
    irdorbs(iorb) = 1
    do ii=1,norb
      orbsao(ii,iorb) = c*rand_cvb(zero)
      if (ii == iorb) orbsao(ii,iorb) = one
    end do
  else
    ! Dummy calls to RAND to get consistent guesses:
    do ii=1,norb
      dum = rand_cvb(zero)
    end do
  end if
end do

! Collect orbitals and transform AO -> MO:
norb_ao = 0
do iorb=1,norb
  if (irdorbs(iorb) == 1) then
    call fmove_cvb(orbsao(1,iorb),orbs(1,iorb),norb)
  else if (irdorbs(iorb) == 2) then
    norb_ao = norb_ao+1
    if (norb_ao /= iorb) call fmove_cvb(orbsao(1,iorb),orbsao(1,norb_ao),nbas_mo)
  end if
end do
i1 = mstackr_cvb(norb*norb_ao)
call ao2mo_cvb(orbsao,work(i1),norb_ao)
iorb_ao = 0
do iorb=1,norb
  if (irdorbs(iorb) == 2) then
    iorb_ao = iorb_ao+1
    call fmove_cvb(work((iorb_ao-1)*norb+i1),orbs(1,iorb),norb)
  end if
end do
call mfreer_cvb(i1)

call nize_cvb(orbs,norb,dum,norb,0,0)

if (abs(detm_cvb(orbs,norb)) < 1d-8) then
  dum = rand_cvb(.777d0)
  c = 1d-1
  do iorb=1,norb
    do ii=1,norb
      orbs(ii,iorb) = orbs(ii,iorb)+c*(one-two*rand_cvb(zero))
    end do
  end do
  if (abs(detm_cvb(orbs,norb)) < 1d-8) then
    if (ip(1) >= 0) write(6,'(a)') ' Starting orbital guess was near-singular - using semi-random guess instead.'
    dum = rand_cvb(.777d0)
    c = 1d-1
    do iorb=1,norb
      do ii=1,norb
        orbs(ii,iorb) = c*rand_cvb(zero)
        if (ii == iorb) orbs(ii,iorb) = one
      end do
    end do
  else
    if (ip(1) >= 0) write(6,'(a)') ' Starting orbital guess was near-singular - scrambling orbital coefficients.'
  end if
  call nize_cvb(orbs,norb,dum,norb,0,0)
end if

call nize_cvb(orbs,norb,dum,norb,0,0)

! Perfect-pairing spin function(s):
if (dnrm2_(nvb,cvb,1) < thresh) then
  kbasiscvb = kbasis
  call ppgs_cvb(cvb)
end if

cnrm = dnrm2_(nvb,cvb,1)
if (cnrm < thresh) then
  write(6,*) ' Fatal error - starting structure coefficients all zero !'
  call abend_cvb()
end if

if (kbasiscvb /= kbasis) then
  call mktrnspn_cvb()
  call untouch_cvb('TRNSPN')
end if

!if(ploc) call rtransf_plc(orbs,cvb)

if ((ip(1) >= 2) .and. (.not. endvar)) then
  write(6,'(/,a)') ' Wavefunction guess :'
  call report_cvb(orbs,norb)
  write(6,'(/,a)') ' Structure coefficients :'
  write(6,'(a)') ' ------------------------'
  call vecprint_cvb(cvb,nvb)
end if

return

end subroutine mkguess2_cvb
