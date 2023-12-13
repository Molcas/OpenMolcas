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

subroutine wrgspr_cvb(recn,c,i,n,ic,ierr)

use casvb_global, only: nbas_mo
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: recn
integer(kind=iwp), intent(in) :: i, n, ic
real(kind=wp), intent(_IN_) :: c(n)
integer(kind=iwp), intent(inout) :: ierr
integer(kind=iwp) :: ioffs, ioffs_cvb, ioffs_orbs, ioffs_orbsao, ioffs_orbslao, kbasiscvb1, nbas_mo1, norb1, nvb1

! Read header:
call rdheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)

if (ic == 1) then
  ! >>> Orbital write, I is orbital number:
  if (i > norb1) then
    ierr = 1
    return
  end if
  ioffs = (i-1)*norb1+ioffs_orbs
  call wrlow_cvb(c,min(norb1,n),recn,ioffs)
else if (ic == 2) then
  ! >>> Structure write, I is starting structure coefficient:
  if (i > nvb1) then
    ierr = 1
    return
  end if
  ioffs = i-1+ioffs_cvb
  call wrlow_cvb(c,min(nvb1,n),recn,ioffs)
else if (ic == 3) then
  ! >>> Write of orbital in AO basis, I is orbital number:
  if (i > norb1) then
    ierr = 1
    return
  end if
  if (nbas_mo1 == 0) then
    nbas_mo1 = nbas_mo
    call wrheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
  end if
  ! Error return if AO bases are not identical:
  if (nbas_mo1 /= n) then
    ierr = 1
    return
  end if
  ioffs = (i-1)*nbas_mo1+ioffs_orbsao
  call wrlow_cvb(c,min(nbas_mo1,n),recn,ioffs)
else if (ic == 4) then
  ! >>> Write of localized orbital in AO basis, I is orbital number:
  if (i > norb1) then
    ierr = 1
    return
  end if
  if (nbas_mo1 == 0) then
    nbas_mo1 = nbas_mo
    call wrheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
  end if
  ! Error return if AO bases are not identical:
  if (nbas_mo1 /= n) then
    ierr = 1
    return
  end if
  ioffs = (i-1)*nbas_mo1+ioffs_orbslao
  call wrlow_cvb(c,min(nbas_mo1,n),recn,ioffs)
end if

return

end subroutine wrgspr_cvb
