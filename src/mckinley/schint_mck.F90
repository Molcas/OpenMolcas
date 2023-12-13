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

subroutine SchInt_mck(CoorM,iAnga,nAlpha,nBeta,nMemab,Zeta,ZInv,rKapab,P,nZeta,Work2,nWork2,Work3,nWork3)
!***********************************************************************
!                                                                      *
! Object: to compute zeta, kappa, P, and the integrals [nm|nm] for     *
!         prescreening. This is done for all unique pairs of centers   *
!         generated from the symmetry unique centers A and B.          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             June '91, modified to compute zeta, P, kappa and inte-   *
!             grals for Schwartz inequality in a k2 loop.              *
!             April '92 modified from k2Loop to a separate subroutine  *
!              for estimates of the gradient.                          *
!***********************************************************************

use Index_Functions, only: nTri3_Elem1, nTri_Elem1
use Real_Spherical, only: ipSph, RSph
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), nAlpha, nBeta, nMemab, nZeta, nWork2, nWork3
real(kind=wp), intent(in) :: CoorM(3,4), Zeta(nZeta), ZInv(nZeta), rKapab(nZeta), P(nZeta,3)
real(kind=wp), intent(inout) :: Work3(nWork3)
real(kind=wp), intent(out) :: Work2(nWork2)
integer(kind=iwp) :: ijklcd, ipIn, la, lb, mabMax, mabMin, mcdMax, mcdMin, mZeta, nijkla, nT
real(kind=wp) :: CoorAC(3,2), Q(3)
logical(kind=iwp) :: NoSpecial
logical(kind=iwp), external :: EQ
external :: TERIS, ModU2, Cff2DS, rys2d

Q(:) = One
la = iAnga(1)
lb = iAnga(2)

! Compute primitive integrals to be used in the prescreening
! by the Schwartz inequality.

! Compute actual size of [a0|c0] block

mabMin = nTri3_Elem1(max(la,lb)-1)
if (EQ(CoorM(1,1),CoorM(1,2))) mabMin = nTri3_Elem1(la+lb-1)
mabMax = nTri3_Elem1(la+lb)-1
mcdMin = mabmin
mcdMax = mabMax

! Find the proper centers to start of with the angular
! momentum on. If la == lb there will exist an
! ambiguity to which center that angular momentum should
! be accumulated on. In that case we will use A and C of
! the order as defined by the basis functions types.

if (iAnga(1) >= iAnga(2)) then
  CoorAC(:,1) = CoorM(:,1)
  CoorAC(:,2) = CoorM(:,3)
else
  CoorAC(:,1) = CoorM(:,2)
  CoorAC(:,2) = CoorM(:,4)
end if

mZeta = nAlpha*nBeta

! Compute [a0|c0], ijkl,a,c

nT = mZeta*1
NoSpecial = .true.
call Rys(iAnga,nT,Zeta,ZInv,mZeta,[One],[One],1,P,nZeta,Q,1,rKapab,[One],CoorM,CoorM,CoorAC,mabMin,mabMax,mcdMin,mcdMax,Work2, &
         nWork2,TERIS,ModU2,Cff2DS,Rys2D,NoSpecial)

! Apply transfer equation to generate [a0|cd], IJKLa,c,d

nijkla = (mabMax-mabMin+1)*mZeta
call HRR(la,lb,CoorM(1,1),CoorM(1,2),Work2,nijkla,nMemab,ipIn)

! Transform to spherical gaussians [a0|CD], CDIJKL,a. This
! will also put the integrals in the right position for the
! transfer equation.

call CrSph1(Work2(ipIn),nijkla,Work3,nWork3,RSph(ipSph(la)),nTri_Elem1(la),nTri_Elem1(la),.false.,RSph(ipSph(lb)),nTri_Elem1(lb), &
            nTri_Elem1(lb),.false.,Work2,nTri_Elem1(la)*nTri_Elem1(lb))

! Apply transfer equation to generate [ab|CD], CDIJKL,a,b

ijklcd = nTri_Elem1(la)*nTri_Elem1(lb)*mZeta
call HRR(la,lb,CoorM(1,1),CoorM(1,2),Work2,ijklcd,nMemab,ipIn)

! Transform to spherical gaussians [AB|CD], IJKL,ABCD

call CrSph2(Work2(ipIn),mZeta,nTri_Elem1(la)*nTri_Elem1(lb),Work3,nWork3,RSph(ipSph(la)),nTri_Elem1(la),nTri_Elem1(la),.false., &
            RSph(ipSph(lb)),nTri_Elem1(lb),nTri_Elem1(lb),.false.,Work2,nTri_Elem1(la)*nTri_Elem1(lb))

return

end subroutine SchInt_mck
