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

subroutine getfree_cvb(nfrr,n_div,nfrdim,iter,fx)

use casvb_global, only: cpu0, cvb, dxmove, formE, icrit, imethod, ipr, nfr, nfrorb, norb, nvb, orbs, proj, projcas, strucopt
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: nfrr, n_div, nfrdim
integer(kind=iwp), intent(in) :: iter
real(kind=wp), intent(in) :: fx
real(kind=wp) :: fxlast = Zero
logical(kind=iwp) :: orb_is_cheap
real(kind=wp), external :: tim_cvb

dxmove = .true.
if (iter >= 0) then
  if (ipr(3) >= 2) then
    write(u6,'(/,a,i5,a,f10.3,a)') ' Iteration',iter,' at',tim_cvb(cpu0),' CPU seconds'
    write(u6,'(a)') ' ---------------------------------------'
  end if
  if (icrit == 1) then
    if (ipr(3) >= 2) write(u6,formE) ' Svb :      ',fx
    if ((ipr(3) >= 2) .and. (iter > 1)) write(u6,formE) ' Svb chg. : ',fx-fxlast
  else if (icrit == 2) then
    if (ipr(3) >= 2) write(u6,formE) ' Evb :      ',fx
    if ((ipr(3) >= 2) .and. (iter > 1)) write(u6,formE) ' Evb chg. : ',fx-fxlast
  end if
  if (ipr(3) >= 2) then
    call report_cvb(orbs,norb)
    if (strucopt) then
      write(u6,'(/,a)') ' Structure coefficients :'
      write(u6,'(a)') ' ------------------------'
      call vecprint_cvb(cvb,nvb)
    end if
  end if
end if
fxlast = fx
call make_cvb('ORBFREE')
call make_cvb('CIFREE')
nfrr = nfr
if (imethod /= 4) then
  nfrdim = max(0,nfr-1)
else
  nfrdim = nfr
end if
orb_is_cheap = ((icrit == 1) .and. (.not. (proj .or. projcas)))
! Set N_DIV:
if ((.not. strucopt) .or. (.not. orb_is_cheap)) then
  n_div = 0
else
  n_div = nfrorb
end if

return

end subroutine getfree_cvb
