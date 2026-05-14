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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PckInt(abab,mZeta,nab,ab,rKappa,Mode,Zeta,nZeta,qKappa)
!***********************************************************************
!                                                                      *
! Object: to keep the diagonal angular indices of a integral batch.    *
!         The integrals are also stripped of the prefactor due to      *
!         the product of gaussians. In case of numerical different-    *
!         iation the prefactor will be due to the undifferentiated     *
!         charge densities.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             April '92                                                *
!***********************************************************************

use Constants, only: Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: mZeta, nab, nZeta
real(kind=wp), intent(in) :: abab(mZeta,nab,nab), rKappa(mZeta), Zeta(mZeta), qKappa(mZeta)
real(kind=wp), intent(inout) :: ab(nZeta,nab)
logical(kind=iwp), intent(in) :: Mode
integer(kind=iwp) :: iab

if (Mode) then
  ! Integrals
  do iab=1,nab
    ab(1:mZeta,iab) = sqrt(sqrt(Two*Zeta(:))*abs(abab(:,iab,iab)))/rKappa(:)
  end do
else
  ! Integrals for numerical estimation of the gradient.
  do iab=1,nab
    ab(1:mZeta,iab) = sqrt(Two*Zeta(:))*abab(:,iab,iab)/(rKappa(:)*qKappa(:))
  end do
end if
#ifdef _DEBUGPRINT_
write(u6,*) 'nZeta,mZeta=',nZeta,mZeta
call RecPrt(' abab','(5G20.10)',abab,mZeta,nab**2)
call RecPrt(' rKappa','(5G20.10)',rKappa,mZeta,1)
call RecPrt(' Zeta  ','(5G20.10)',Zeta,mZeta,1)
do iab=1,nab
  call RecPrt(' ab ','(5G20.10)',ab(1,iab),mZeta,1)
end do
#endif

return

end subroutine PckInt
