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
! Copyright (C) 1990-1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DoZeta(Alpha,nAlpha,Beta,nBeta,A,B,P,Zeta,rKappa,ZInv,Alpha_,Beta_,Ind_Pair)
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

use Constants, only: One, TwoP54
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAlpha, nBeta
real(kind=wp), intent(in) :: Alpha(nAlpha), Beta(nBeta), A(3), B(3)
real(kind=wp), intent(out) :: P(nAlpha*nBeta,3), Zeta(nAlpha*nBeta), rKappa(nAlpha*nBeta), ZInv(nAlpha*nBeta), &
                              Alpha_(nAlpha*nBeta), Beta_(nAlpha*nBeta)
integer(kind=iwp), intent(out) :: Ind_Pair(nAlpha*nBeta+1)
integer(kind=iwp) :: iAlpha, iBeta, iZeta
real(kind=wp) :: AB2, Tmp0, Tmp1

#ifdef _DEBUGPRINT_
call RecPrt(' In DoZeta:Alpha',' ',Alpha,nAlpha,1)
call RecPrt(' In DoZeta:Beta',' ',Beta,nBeta,1)
#endif

AB2 = (A(1)-B(1))**2+(A(2)-B(2))**2+(A(3)-B(3))**2
do iBeta=1,nBeta
  do iAlpha=1,nAlpha
    iZeta = nAlpha*(iBeta-1)+iAlpha
    Zeta(iZeta) = Alpha(iAlpha)+Beta(iBeta)
    Alpha_(iZeta) = Alpha(iAlpha)
    Beta_(iZeta) = Beta(iBeta)
    ZInv(iZeta) = One/Zeta(iZeta)
    Tmp0 = ZInv(iZeta)
    Tmp1 = TwoP54*exp(-Alpha(iAlpha)*Beta(iBeta)*AB2*Tmp0)*Tmp0
    if (Tmp1 < 1.0e-99_wp) Tmp1 = 1.0e-99_wp
    rKappa(iZeta) = Tmp1
    P(iZeta,1) = (Alpha(iAlpha)*A(1)+Beta(iBeta)*B(1))*Tmp0
    P(iZeta,2) = (Alpha(iAlpha)*A(2)+Beta(iBeta)*B(2))*Tmp0
    P(iZeta,3) = (Alpha(iAlpha)*A(3)+Beta(iBeta)*B(3))*Tmp0
    Ind_Pair(iZeta) = iZeta
  end do
end do
Ind_Pair(nAlpha*nBeta+1) = nAlpha*nBeta

! Sort from Large to Small

!#define _New_Code_
#ifdef _New_Code_
do iZeta=1,nAlpha*nBeta-1
  Tmp1 = rKappa(iZeta)
  do jZeta=iZeta+1,nAlpha*nBeta
    if (Tmp1 < rKappa(jZeta)) then
      Tmp1 = rKappa(jZeta)
      Tmp2 = Zeta(iZeta)
      Zeta(iZeta) = Zeta(jZeta)
      Zeta(jZeta) = Tmp2
      Tmp2 = Alpha_(iZeta)
      Alpha_(iZeta) = Alpha_(jZeta)
      Alpha_(jZeta) = Tmp2
      Tmp2 = Beta_(iZeta)
      Beta_(iZeta) = Beta_(jZeta)
      Beta_(jZeta) = Tmp2
      Tmp2 = ZInv(iZeta)
      ZInv(iZeta) = ZInv(jZeta)
      ZInv(jZeta) = Tmp2
      Tmp2 = rKappa(iZeta)
      rKappa(iZeta) = rKappa(jZeta)
      rKappa(jZeta) = Tmp2
      Tmp2 = P(iZeta,1)
      P(iZeta,1) = P(jZeta,1)
      P(jZeta,1) = Tmp2
      Tmp2 = P(iZeta,2)
      P(iZeta,2) = P(jZeta,2)
      P(jZeta,2) = Tmp2
      Tmp2 = P(iZeta,3)
      P(iZeta,3) = P(jZeta,3)
      P(jZeta,3) = Tmp2
      iTmp = Ind_Pair(iZeta)
      Ind_Pair(iZeta) = Ind_Pair(jZeta)
      Ind_Pair(jZeta) = iTmp
    end if
  end do
end do
#endif

#ifdef _DEBUGPRINT_
call RecPrt(' In DoZeta: Kappa',' ',rKappa,nAlpha*nBeta,1)
call RecPrt(' In DoZeta: P',' ',P,nAlpha*nBeta,3)
#endif

return

end subroutine DoZeta
