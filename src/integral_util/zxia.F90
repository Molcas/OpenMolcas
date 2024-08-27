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
! Copyright (C) 1990,2020, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ZXia(Zeta,ZInv,N,M,Alpha,Beta)
!***********************************************************************
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January 1990                                             *
!***********************************************************************

use Constants, only: One

implicit none
integer, intent(In) :: N, M
real*8, intent(In) :: Alpha(N), Beta(M)
real*8, intent(InOut) :: Zeta(N,M), ZInv(N,M)
integer j

#ifdef _DEBUGPRINT_
call RecPrt(' In ZXia: Alpha',' ',Alpha,N,1)
call RecPrt(' In ZXia: Beta ',' ',Beta,M,1)
#endif

do j=1,M
  Zeta(:,j) = (Alpha(:)+Beta(j))
  ZInv(:,j) = One/Zeta(:,j)
end do
#ifdef _DEBUGPRINT_
call RecPrt(' In ZXia: Zeta',' ',Zeta,N,M)
#endif

return

end subroutine ZXia
