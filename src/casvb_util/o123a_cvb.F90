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

subroutine o123a_cvb( &
#                    define _CALLING_
#                    include "opta_interface.fh"
                    )

use casvb_global, only: eigval, eigvec, ip, ograd, ogradp
use Definitions, only: iwp, u6

implicit none
#include "opta_interface.fh"

call gethess_cvb(eigvec)
call mxdiag_cvb(eigvec,eigval,nparam)
call mxatb_cvb(ograd,eigvec,1,nparam,nparam,ogradp)
if (ip >= 2) then
  write(u6,'(a)') ' Gradient in basis of Hessian eigenvectors :'
  call vecprint_cvb(ogradp,nparam)
end if

return

end subroutine o123a_cvb
