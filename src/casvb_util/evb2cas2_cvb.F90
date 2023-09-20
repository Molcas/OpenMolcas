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

subroutine evb2cas2_cvb(orbs,cvb,ioptc,iter,fx,dxnrm,dx_amx,civec,civb,civbh,res,resh,cvbdet,gjorb)

use casvb_global, only: dx, formAD, formAF, grd

implicit real*8(a-h,o-z)
! ... Files/Hamiltonian available ...
logical, external :: tstfile_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
logical dx_ok, grad_ok
dimension orbs(norb,norb), cvb(nvb)
dimension civec(ndet), civb(ndet), civbh(ndet)
dimension cvbdet(ndetvb)
dimension gjorb(*), res(*), resh(*)
dimension h(2,2), eig(2)
dimension orbinv(norb,norb)

if (ip(3) >= 0) then
  write(6,'(/,a)') ' Starting VB2CAS optimization.'
  write(6,'(a)') ' -----------------------------'
end if

dx_ok = ((dx_amx < dx(1,3)) .and. (dxnrm < dx(2,3)))

if (.not. projcas) then
  call str2vbc_cvb(cvb,cvbdet)
  call vb2cic_cvb(cvbdet,civb)
else if (projcas) then
  if (memplenty) then
    call cicopy_cvb(civec,civbh)
  else
    call cird_cvb(civbh,61001.2d0)
  end if
  call fmove_cvb(orbs,orbinv,norb*norb)
  call mxinv_cvb(orbinv,norb)
  call gaussj_cvb(orbinv,gjorb)
  call applyt_cvb(civbh,gjorb)
  call pvbcopy_cvb(civbh,civb)
  call ci2vbc_cvb(civbh,cvbdet)
end if

call gaussj_cvb(orbs,gjorb)
call applyt_cvb(civb,gjorb)
call proj_cvb(civb)

call cinorm_cvb(civb,cnrm)
call ciscale_cvb(civb,one/sqrt(cnrm))

call cicopy_cvb(civb,civbh)
call applyh_cvb(civbh)

call cidot_cvb(civb,civbh,evb)
if (ip(3) >= 2) write(6,formAF) ' Residual calculation based on Evb :',evb+corenrg
! RES()=CIVBH()-EVB*CIVB()
call cicopy_cvb(civbh,res)
call cidaxpy_cvb(-evb,civb,res)

if (tstfile_cvb(67123.2d0)) then
  call cird_cvb(resh,67123.2d0)
  call cidot_cvb(res,resh,rescas_ovr)
  !call cidot_cvb(civb,resh,civb_ovr)
  !dxnrm_ci = sqrt(2d0*(one-civb_ovr))
  !write(6,*) ' dxnrm dxnrm_ci :',dxnrm,dxnrm_ci
  !write(6,*) ' gradient in VB basis :',2d0*rescas_ovr/dxnrm
  grad_ok = (2d0*rescas_ovr/dxnrm < grd(1,3))
else
  grad_ok = .false.
end if
call ciwr_cvb(civb,67123.2d0)

call cinorm_cvb(res,resnrm)
if (ip(3) >= 2) then
  write(6,'(a)') ' '
  write(6,formAD) ' Residual norm:',resnrm
  write(6,'(a)') ' '
end if
call ciscale_cvb(res,one/sqrt(resnrm))
call cidot_cvb(res,civb,ovr)
! RES()=RES()-OVR*CIVB()
call cidaxpy_cvb(-ovr,civb,res)
call cinorm_cvb(res,resnrm)
call ciscale_cvb(res,one/sqrt(resnrm))
call cidot_cvb(civbh,civb,h(1,1))
call cidot_cvb(civbh,res,h(1,2))

call cicopy_cvb(res,resh)
call applyh_cvb(resh)

call cidot_cvb(resh,civb,h(2,1))
call cidot_cvb(resh,res,h(2,2))

if (ip(3) >= 2) then
  write(6,*) ' 2x2 Hamiltonian matrix :'
  eig(1) = h(1,1)
  eig(2) = h(2,2)
  h(1,1) = h(1,1)+corenrg
  h(2,2) = h(2,2)+corenrg
  call mxprintd_cvb(h,2,2,0)
  h(1,1) = eig(1)
  h(2,2) = eig(2)
end if

call mxdiag_cvb(h,eig,2)
if (ip(3) >= 2) then
  write(6,*) ' Eigenvalues :',eig(1)+corenrg,eig(2)+corenrg
  write(6,*) ' Eigenvectors :'
  call mxprint_cvb(h,2,2,0)
end if

if (abs(h(1,1)) > abs(h(1,2))) then
  if (ip(3) >= 2) write(6,*) ' Using root 1 :'
  call ciscale_cvb(civb,h(1,1))
  call cidaxpy_cvb(h(2,1),res,civb)
else
  if (ip(3) >= 2) write(6,*) ' Using root 2 :'
  call ciscale_cvb(civb,h(1,2))
  call cidaxpy_cvb(h(2,2),res,civb)
end if
call cinorm_cvb(civb,cnrm)
call ciscale_cvb(civb,one/sqrt(cnrm))
if (memplenty) then
  !call cidot_cvb(civb,civec,ovr)
  call cicopy_cvb(civb,civec)
else
  call cird_cvb(res,61001.2d0)
  !call cidot_cvb(civb,res,ovr)
  call ciwr_cvb(civb,61001.2d0)
end if

evb = evb+corenrg
fx = evb
ovraa = one
iter = 1
ioptc = 0
!ovrcrit = .125d-10
!if (abs(one-abs(ovr)) > ovrcrit) iter = 2
if (.not. (dx_ok .and. grad_ok)) iter = 2
call setcnt_cvb(civec,1)

return

end subroutine evb2cas2_cvb
