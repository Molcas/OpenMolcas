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

subroutine testconv2_cvb(close2conv,converged,wrongstat,act,zz,step,grad,npr,eigmn,eigmx,eigmna,nposeig,nnegeig)

use casvb_global, only: dx, dfx, formAD, grd, ipr, sgn, singul, zzmin, zzmax
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(inout) :: close2conv
logical(kind=iwp), intent(out) :: converged, wrongstat
integer(kind=iwp), intent(in) :: npr, nposeig, nnegeig
real(kind=wp), intent(in) :: act, zz, step(npr), grad(npr), eigmn, eigmx, eigmna
integer(kind=iwp) :: idum, mm
real(kind=wp) :: grad_amx, grad_nrm, grad_rms, step_amx, step_nrm, step_rms
logical(kind=iwp) :: close_old, correct_index, dfx_is_small, grad_is_small, step_is_small, zz_ok
real(kind=wp), external :: dnrm2_

close_old = close2conv

step_nrm = dnrm2_(npr,step,1)
step_rms = step_nrm/sqrt(real(npr,kind=wp))
call findamx_cvb(step,npr,step_amx,idum)

grad_nrm = dnrm2_(npr,grad,1)
grad_rms = grad_nrm/sqrt(real(npr,kind=wp))
call findamx_cvb(grad,npr,grad_amx,idum)

if (ipr(3) >= 2) then
  if ((nnegeig > 0) .and. (nposeig > 0)) then
    write(u6,*) ' Maximum eigenvalue : ',eigmx,' of first ',nnegeig,' values.'
    write(u6,*) ' Minimum eigenvalue : ',eigmn,' of last  ',nposeig,' values.'
  else if (nnegeig > 0) then
    write(u6,*) ' Maximum eigenvalue : ',eigmx
  else
    write(u6,*) ' Minimum eigenvalue : ',eigmn
  end if
end if

! Local region / final convergence reached if:
! a) Change in F(X) small
! b) Step is small
! c) Gradient is small
! d) Hessian has correct index
! e) ACT/EXP ratio is within acceptable interval

if (eigmna > singul(1)) then
  mm = 1
else
  mm = 2
end if
dfx_is_small = act < dfx(mm)
step_is_small = (step_nrm < dx(2,mm)) .and. (step_rms < dx(3,mm)) .and. (step_amx < dx(1,mm))
grad_is_small = (grad_nrm < grd(2,mm)) .and. (grad_rms < grd(3,mm)) .and. (grad_amx < grd(1,mm))
correct_index = (eigmx < sgn(mm)) .and. (eigmn > -sgn(mm))
zz_ok = (zz > zzmin(mm)) .and. (zz < zzmax(mm))

close2conv = dfx_is_small .and. step_is_small .and. grad_is_small .and. correct_index .and. zz_ok

if (eigmna > singul(2)) then
  mm = 3
else
  mm = 4
end if
dfx_is_small = act < dfx(mm)
step_is_small = (step_nrm < dx(2,mm)) .and. (step_rms < dx(3,mm)) .and. (step_amx < dx(1,mm))
grad_is_small = (grad_nrm < grd(2,mm)) .and. (grad_rms < grd(3,mm)) .and. (grad_amx < grd(1,mm))
correct_index = (eigmx < sgn(mm)) .and. (eigmn > -sgn(mm))
zz_ok = (zz > zzmin(mm)) .and. (zz < zzmax(mm))

if (ipr(3) >= 2) then
  write(u6,'(/,a)') ' Test of convergence :'
  write(u6,'(a)') ' ---------------------'
  call cvprt_cvb(' 1) Change in F(x) :',dfx_is_small)
  call cvprt_cvb(' 2) Step length    :',step_is_small)
  call cvprt_cvb(' 3) Grad norm      :',grad_is_small)
  call cvprt_cvb(' 4) Hessian index  :',correct_index)
  call cvprt_cvb(' 5) Act/Exp ratio  :',zz_ok)
  write(u6,*) ' '
  call cvprt2_cvb(' F(x) change   :',act,dfx(mm),1)
  call cvprt2_cvb(' Norm of step  :',step_nrm,dx(2,mm),1)
  call cvprt2_cvb(' RMS of step   :',step_rms,dx(3,mm),1)
  call cvprt2_cvb(' AMAX of step  :',step_amx,dx(1,mm),1)
  call cvprt2_cvb(' Norm of grad  :',grad_nrm,grd(2,mm),1)
  call cvprt2_cvb(' RMS of grad   :',grad_rms,grd(3,mm),1)
  call cvprt2_cvb(' AMAX of grad  :',grad_amx,grd(1,mm),1)
  call cvprt2_cvb(' Max. eigval   :',eigmx,sgn(mm),1)
  call cvprt2_cvb(' Min. eigval   :',eigmn,-sgn(mm),2)
  call cvprt2_cvb(' Act/Exp ratio :',zz,zzmin(mm),2)
  call cvprt2_cvb(' Act/Exp ratio :',zz,zzmax(mm),1)
  write(u6,*) ' '
end if

converged = dfx_is_small .and. step_is_small .and. grad_is_small .and. correct_index .and. zz_ok
converged = converged .and. close2conv

if (eigmna > singul(3)) then
  mm = 5
else
  mm = 6
end if
dfx_is_small = act < dfx(mm)
step_is_small = (step_nrm < dx(2,mm)) .and. (step_rms < dx(3,mm)) .and. (step_amx < dx(1,mm))
grad_is_small = (grad_nrm < grd(2,mm)) .and. (grad_rms < grd(3,mm)) .and. (grad_amx < grd(1,mm))
correct_index = (eigmx < sgn(mm)) .and. (eigmn > -sgn(mm))
zz_ok = (zz > zzmin(mm)) .and. (zz < zzmax(mm))

! Wrong stationary point if otherwise converged, but Hessian
! hasn't got correct index:
wrongstat = dfx_is_small .and. step_is_small .and. grad_is_small .and. (.not. correct_index) .and. zz_ok

if ((ipr(3) >= 1) .and. close2conv .and. (.not. (converged .or. close_old))) write(u6,'(a)') ' Optimization entering local region.'
if (converged .and. (ipr(3) >= 0)) then
  write(u6,formAD) ' Converged ... maximum update to coefficient:',step_amx
  if (eigmna <= singul(2)) then
    write(u6,'(a)') ' Warning - singular hessian!'
    write(u6,formAD) ' Smallest Hessian eigenvalue :',eigmna
  end if
end if

return

end subroutine testconv2_cvb
