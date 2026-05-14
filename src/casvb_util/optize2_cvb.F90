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

subroutine optize2_cvb(fx,nparm,ioptc,iter_is_1,opta,optb)

use casvb_global, only: endwhenclose, expct, formAD, formAF, fxbest, hh, ip, maxize, odx, ograd
use casvb_interfaces, only: opta_sub, optb_sub
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: fx
integer(kind=iwp), intent(in) :: nparm
integer(kind=iwp), intent(inout) :: ioptc
logical(kind=iwp), intent(in) :: iter_is_1
procedure(opta_sub) :: opta
procedure(optb_sub) :: optb
integer(kind=iwp) :: iopth
real(kind=wp) :: dxnrm, exp_tc, grdnrm, s11, s12, s22
logical(kind=iwp) :: close2conv = .false., close2conv_begin, converged, first_time, opth, scalesmall1, skipupd, wrongstat
real(kind=wp), external :: ddot_, dnrm2_

converged = .false.
if (iter_is_1) close2conv = .false.

! << Initialization >>
call grad_cvb(ograd)
call ddproj_cvb(ograd,nparm)
grdnrm = dnrm2_(nparm,ograd,1)
call opta(nparm)

if (ip >= 2) write(u6,'(/a)') ' *****   2. order optimizer   *****'

! << Now trust region control >>
exp_tc = expct
first_time = .true.
opth = .false.
iopth = 0
do
  call trust_cvb(iopth,opth,maxize,fx,fxbest,expct,hh,dxnrm,ioptc,scalesmall1,close2conv,converged,skipupd)
  if (ioptc == -2) return

  ! << Make update >>
  if (.not. (skipupd .or. (hh == Zero))) then

    do
      close2conv_begin = close2conv

      call optb(nparm,dxnrm,grdnrm,close2conv)

      if (.not. first_time) exit
      first_time = .false.

      call testconv_cvb(fx,nparm,odx,ograd,exp_tc,close2conv,converged,wrongstat)

      if ((.not. close2conv) .or. close2conv_begin) exit
    end do
    if (((ip == 2) .and. (.not. opth)) .or. (ip >= 3)) then
      s11 = ddot_(nparm,odx,1,odx,1)
      s22 = ddot_(nparm,ograd,1,ograd,1)
      s12 = ddot_(nparm,odx,1,ograd,1)
      write(u6,formAD) ' Overlap between normalized vectors <DX|GRAD> :',s12/sqrt(s11*s22)
    end if
    call fxdx_cvb(fx,.false.,odx)
  end if
  if (.not. opth) exit
end do

if ((ioptc >= -1) .and. (hh /= Zero)) then
  if (ip >= 2) then
    write(u6,'(a)') ' '
    write(u6,formAF) ' HH & norm of update :',hh,dxnrm
  end if
  call update_cvb(odx)
end if
if (converged) then
  ioptc = 0
else if (close2conv .and. endwhenclose) then
  ioptc = -3
else
  ioptc = 1
end if

return

end subroutine optize2_cvb
