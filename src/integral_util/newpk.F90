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
! Copyright (C) 1990-1992,1999, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NewPK(A,B,P,mZeta,nZeta,rKappa,Alpha,Beta)
!***********************************************************************
!                                                                      *
! Object : to compute P and kappa.                                     *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!             May '90, modified for integral cutoff.                   *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN.                              *
!             June '91, modified for k2 loop.                          *
!             January '92, modified for gradient calculations.         *
!***********************************************************************

use Constants, only: Zero, One, TwoP54
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mZeta, nZeta
real(kind=wp), intent(in) :: A(3), B(3), Alpha(nZeta), Beta(nZeta)
real(kind=wp), intent(out) :: P(nZeta,3), rKappa(nZeta)
integer(kind=iwp) :: iZeta
real(kind=wp) :: AB2, Tmp0, Tmp1

#ifdef _DEBUGPRINT_
call RecPrt(' In NewPK:Alpha',' ',Alpha,mZeta,1)
call RecPrt(' In NewPK:Beta',' ',Beta,mZeta,1)
#endif
AB2 = (A(1)-B(1))**2+(A(2)-B(2))**2+(A(3)-B(3))**2
do iZeta=1,mZeta
  Tmp0 = One/(Alpha(iZeta)+Beta(iZeta))
  Tmp1 = TwoP54*exp(-Alpha(iZeta)*Beta(iZeta)*AB2*Tmp0)*Tmp0
  if (Tmp1 < 1.0e-99_wp) Tmp1 = 1.0e-99_wp
  rKappa(iZeta) = Tmp1
  P(iZeta,:) = (Alpha(iZeta)*A(:)+Beta(iZeta)*B(:))*Tmp0
end do
rKappa(mZeta+1:) = Zero
P(mZeta+1:,:) = Zero

#ifdef _DEBUGPRINT_
call RecPrt(' In NewPK: Kappa',' ',rKappa,mZeta,1)
call RecPrt(' In NewPK: Px',' ',P(1,1),mZeta,1)
call RecPrt(' In NewPK: Py',' ',P(1,2),mZeta,1)
call RecPrt(' In NewPK: Px',' ',P(1,3),mZeta,1)
#endif

return

end subroutine NewPK
