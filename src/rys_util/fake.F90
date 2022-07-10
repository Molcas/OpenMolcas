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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine Fake(U2,mT,nRys,ZEInv)
!***********************************************************************
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             May '90                                                  *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mT, nRys
real(kind=wp), intent(inout) :: U2(mT,nRys)
real(kind=wp), intent(in) :: ZEInv(mT)

#include "macros.fh"
unused_var(U2)
unused_var(ZEInv)

return

end subroutine Fake
