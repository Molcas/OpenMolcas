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

subroutine o8b_cvb( &
#                  define _CALLING_
#                  include "optb_interface.fh"
                  )

use casvb_global, only: eigval, eigvec, hh, ip, odx, ograd, scalesmall
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "optb_interface.fh"
integer(kind=iwp) :: iprm, ipu, iroot
real(kind=wp) :: fac1
real(kind=wp), external :: dnrm2_

#include "macros.fh"
unused_var(grdnrm)

eigvec(1:nparm+1,1:nparm+1) = Zero
eigvec(2:nparm+1,1) = ograd(1:nparm)
eigvec(1,2:nparm+1) = ograd(1:nparm)
do iprm=2,nparm+1
  eigvec(iprm,iprm) = One
  call hess_cvb(eigvec(2,iprm))
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
odx(:) = eigvec(2:,iroot)
if (abs(eigvec(1,iroot)) > 1.0e-8_wp) then
  fac1 = One/eigvec(1,iroot)
else
  fac1 = sign(One,eigvec(1,iroot))
end if
odx(1:nparm) = fac1*odx(1:nparm)
dxnrm = dnrm2_(nparm,odx,1)
if (.not. close2conv) then
  ipu = 1
else
  ipu = 2
end if
if ((dxnrm > hh) .or. scalesmall(ipu)) then
  odx(1:nparm) = hh/dxnrm*odx(1:nparm)
  dxnrm = hh
end if

return

end subroutine o8b_cvb
