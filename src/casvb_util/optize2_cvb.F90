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

subroutine optize2_cvb(fx,nparm,ioptc,dx,grad,iter_is_1,opta,optb)

use casvb_global, only: endwhenclose, expct, formAD, formAF, fxbest, hh, ip, maxize

implicit real*8(a-h,o-z)
logical opth, skipupd, first_time
logical iter_is_1, close2conv_begin
logical close2conv, converged, wrongstat, scalesmall1
external opta, optb
dimension dx(nparm), grad(nparm)
save zero, close2conv, converged
data zero/0.d0/

converged = .false.
if (iter_is_1) close2conv = .false.

! << Initialization >>
call grad_cvb(grad)
call ddproj_cvb(grad,nparm)
grdnrm = dnrm2_(nparm,grad,1)
call opta(nparm)

if (ip >= 2) write(6,'(/a)') ' *****   2. order optimizer   *****'

! << Now trust region control >>
exp_tc = expct
first_time = .true.
opth = .false.
iopth = 0
do
  call trust_cvb(iopth,opth,maxize,fx,fxbest,expct,hh,dxnrm,ioptc,scalesmall1,close2conv,converged,skipupd)
  if (ioptc == -2) return

  ! << Make update >>
  if (.not. (skipupd .or. (hh == zero))) then

    do
      close2conv_begin = close2conv

      call optb(nparm,dxnrm,grdnrm,close2conv)

      if (.not. first_time) exit
      first_time = .false.

      call testconv_cvb(fx,nparm,dx,grad,exp_tc,close2conv,converged,wrongstat)

      if ((.not. close2conv) .or. close2conv_begin) exit
    end do
    if (((ip == 2) .and. (.not. opth)) .or. (ip >= 3)) then
      s11 = ddot_(nparm,dx,1,dx,1)
      s22 = ddot_(nparm,grad,1,grad,1)
      s12 = ddot_(nparm,dx,1,grad,1)
      write(6,formAD) ' Overlap between normalized vectors <DX|GRAD> :',s12/sqrt(s11*s22)
    end if
    call fxdx_cvb(fx,.false.,dx)
  end if
  if (.not. opth) exit
end do

if ((ioptc >= -1) .and. (hh /= zero)) then
  if (ip >= 2) then
    write(6,'(a)') ' '
    write(6,formAF) ' HH & norm of update :',hh,dxnrm
  end if
  call update_cvb(dx)
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
