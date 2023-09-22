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

subroutine opt2_cvb(orbs,cvb)

use casvb_global, only: formE
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
real(kind=wp) :: orbs(norb,norb), cvb(nvb)
#include "optze_cvb.fh"
#include "print_cvb.fh"
integer(kind=iwp) :: ioptc, iter
real(kind=wp) :: fx

! The following initialization is to appease a compiler
fx = Zero
ioptc = 0
iter = 0

if (imethod == 11) then
  ! Method = None:
else if (imethod == 4) then
  if (icrit == 1) then
    call svbd_cvb(orbs,cvb,fx,ioptc,iter)
  else if (icrit == 2) then
    call evbd_cvb(orbs,cvb,fx,ioptc,iter)
  end if
else if (imethod == 6) then
  call evb2cas_cvb(orbs,cvb,fx,ioptc,iter)
else

  fx = Zero

  call optize_cvb(fx,ioptc,iter,imethod,isaddle,mxiter,icrit == 1,corenrg,ip(3),ip(4)-2,ip(4)-2,strucopt)

  if ((ioptc == -1) .and. (mxiter > 0)) then
    if (ip(3) >= 0) write(u6,'(a,i4)') ' Maximum number of iterations reached:',mxiter
    if (ip(3) >= 0) write(u6,'(a)') ' Calculation NOT converged!!!'
  end if
end if
if (icrit == 1) then
  svb = fx
else
  evb = fx
end if
if (ip(5) >= 0) then
  if (imethod /= 11) then
    if (icrit == 1) write(u6,formE) ' Final Svb :',svb
    if (icrit == 2) write(u6,formE) ' Final Evb :',evb
  end if
  if ((ip(3) <= 1) .and. (ioptc /= -1)) write(u6,'(a,i4)') ' Number of iterations used:',iter
end if
if (ip(5) >= 2) then
  call report_cvb(orbs,norb)
  write(u6,'(/,a)') ' Structure coefficients :'
  write(u6,'(a)') ' ------------------------'
  call vecprint_cvb(cvb,nvb)
end if
convinone = ((ioptc == 0) .and. (iter <= 1)) .or. endvar
ioptc_new = ioptc
n_iter = n_iter+iter
if (ioptc == 1) ioptc_new = mxiter
if (ioptc == 0) ioptc_new = iter

return

end subroutine opt2_cvb
