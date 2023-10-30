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

subroutine mksymorbs_cvb()

use casvb_global, only: ipr, nconstr, norb, orbs, sorbs, sym
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: ierr, nconstr_kp
real(kind=wp) :: delorbs, dum(1)
real(kind=wp), parameter :: thresh = 1.0e-7_wp
real(kind=wp), external :: detm_cvb, dnrm2_

if (sym) then
  sorbs(:,:) = orbs(:,:)
  nconstr_kp = nconstr
  nconstr = 0
  call symtrizorbs_cvb(orbs)
  nconstr = nconstr_kp
  sorbs(:,:) = orbs(:,:)-sorbs(:,:)
  delorbs = dnrm2_(norb*norb,sorbs,1)
  if ((delorbs > thresh) .and. (ipr(1) >= 2)) then
    write(u6,'(/,a)') ' Change in symmetrized orbitals:'
    call report_cvb(sorbs,norb)
  end if
  ierr = 0
  call nize_cvb(orbs,norb,dum,norb,0,ierr)
  if ((delorbs > thresh) .and. (ipr(1) >= 2)) then
    write(u6,'(a)') ' Orbitals after symmetrization:'
    call report_cvb(orbs,norb)
  end if
  if (abs(detm_cvb(orbs,norb)) < 1.0e-8_wp) then
    write(u6,*) ' Fatal error - orbital matrix singular after symmetrization!'
    call abend_cvb()
  end if
end if

return

end subroutine mksymorbs_cvb
