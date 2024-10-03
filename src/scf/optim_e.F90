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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************

function optim_E(C,G,H,n)

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: optim_E
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: C(n), G(n), H(n,n)
integer(kind=iwp) :: k

Optim_E = sum(C(:)*G(:))
do k=1,n
  Optim_E = Optim_E+Half*C(k)*sum(C(:)*H(k,:))
end do

end function optim_E
