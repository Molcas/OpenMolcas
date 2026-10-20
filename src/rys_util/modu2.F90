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

subroutine ModU2( &
#                define _CALLING_
#                include "modu2_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: precompute u2/(zeta+eta)                                     *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             May '90                                                  *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
#include "modu2_interface.fh"
integer(kind=iwp) :: iT

#ifdef _DEBUGPRINT_
call RecPrt(' In ModU2: U2',' ',U2,nRys,mT)
call RecPrt(' In ModU2: ZEInv',' ',ZEInv,1,mT)
#endif

if (nRys > 1) then
  do iT=1,mT
    U2(:,iT) = U2(:,iT)*ZEInv(iT)
  end do
else
  U2(1,:) = U2(1,:)*ZEInv(:)
end if

end subroutine ModU2
