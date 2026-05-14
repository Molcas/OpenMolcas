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

subroutine rdgspr_cvb(recn,c,i,n,ic,ierr)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: recn
integer(kind=iwp), intent(in) :: i, n, ic
real(kind=wp), intent(out) :: c(n)
integer(kind=iwp), intent(out) :: ierr
integer(kind=iwp) :: ioffs, ioffs_cvb, ioffs_orbs, ioffs_orbsao, ioffs_orbslao, kbasiscvb1, nbas_mo1, norb1, nvb1

ierr = 0
c(:) = Zero
! Read header:
call rdheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)

if (ic == 1) then
  ! >>> Orbital read, I is orbital number:
  if (i > norb1) then
    ierr = 1
    return
  end if
  ioffs = (i-1)*norb1+ioffs_orbs
  call rdlow_cvb(c,min(norb1,n),recn,ioffs)
else if (ic == 2) then
  ! >>> Structure read, I is starting structure coefficient:
  if (i > nvb1) then
    ierr = 1
    return
  end if
  ioffs = i-1+ioffs_cvb
  call rdlow_cvb(c,min(nvb1,n),recn,ioffs)
else if (ic == 3) then
  ! >>> Read of orbital in AO basis, I is orbital number:
  if (i > norb1) then
    ierr = 1
    return
  end if
  ! Error return if AO bases are not identical:
  if (nbas_mo1 /= n) then
    ierr = 1
    return
  end if
  ioffs = (i-1)*nbas_mo1+ioffs_orbsao
  call rdlow_cvb(c,min(nbas_mo1,n),recn,ioffs)
else if (ic == 4) then
  ! >>> Read of localized orbital in AO basis, I is orbital number:
  if (i > norb1) then
    ierr = 1
    return
  end if
  ! Error return if AO bases are not identical:
  if (nbas_mo1 /= n) then
    ierr = 1
    return
  end if
  ioffs = (i-1)*nbas_mo1+ioffs_orbslao
  call rdlow_cvb(c,min(nbas_mo1,n),recn,ioffs)
end if

return

end subroutine rdgspr_cvb
