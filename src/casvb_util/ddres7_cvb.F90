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

subroutine ddres7_cvb( &
#                     define _CALLING_
#                     include "ddres_interface.fh"
                     )

use casvb_global, only: iroot, jroot
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "ddres_interface.fh"
integer(kind=iwp) :: i

#include "macros.fh"
unused_var(rhs)

res(:) = Zero
do i=1,itdav
  res(:) = res(:)+(axc(:,i)-eig_res*sxc(:,i))*solp_res(i)
end do

is_converged = (jroot == iroot)

return

end subroutine ddres7_cvb
