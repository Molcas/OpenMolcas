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

subroutine o8b2_cvb(nparm,dx,grad,eigvec,eigval,dxnrm,close2conv)

use casvb_global, only: hh, ip, scalesmall
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nparm
real(kind=wp) :: dx(nparm), grad(nparm), eigvec(nparm+1,nparm+1), eigval(nparm+1), dxnrm
logical(kind=iwp) :: close2conv
integer(kind=iwp) :: iprm, ipu, iroot
real(kind=wp) :: fac1
real(kind=wp), external :: dnrm2_

call fzero(eigvec,(nparm+1)*(nparm+1))
do iprm=1,nparm
  eigvec(iprm+1,1) = grad(iprm)
  eigvec(1,iprm+1) = grad(iprm)
  eigvec(iprm+1,iprm+1) = One
  call hess_cvb(eigvec(2,iprm+1))
end do
write(u6,*) ' Augmented Hessian matrix :'
call mxprint_cvb(eigvec,nparm+1,nparm+1,0)
call mxdiag_cvb(eigvec,eigval,nparm+1)
iroot = nparm+1
if (ip >= 2) then
  write(u6,'(a)') ' Eigenvalues of augmented Hessian :'
  call vecprint_cvb(eigval,nparm+1)
  write(u6,'(a)') ' Eigenvector to be followed :'
  call vecprint_cvb(eigvec(1,iroot),nparm+1)
end if
write(u6,*) ' Following root no :',iroot
call fmove_cvb(eigvec(2,iroot),dx,nparm)
if (abs(eigvec(1,iroot)) > 1.0e-8_wp) then
  fac1 = One/eigvec(1,iroot)
else
  fac1 = sign(One,eigvec(1,iroot))
end if
call dscal_(nparm,fac1,dx,1)
dxnrm = dnrm2_(nparm,dx,1)
if (.not. close2conv) then
  ipu = 1
else
  ipu = 2
end if
if ((dxnrm > hh) .or. scalesmall(ipu)) then
  call dscal_(nparm,hh/dxnrm,dx,1)
  dxnrm = hh
end if

return

end subroutine o8b2_cvb
