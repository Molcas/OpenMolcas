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
subroutine SchInt(CoorM,iAnga,mZeta,Zeta,ZInv,rKapab,P,rKapcd,Q,nZeta,Wrk,nWork2,HMtrx,nHrrMtrx,iShlla,jShllb,i_Int)
!***********************************************************************
!                                                                      *
! Object: to compute zeta, kappa, P, and the integrals [nm|nm] for     *
!         prescreening. This is done for all unique pairs of centers   *
!         generated from the symmetry unique centers A and B.          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN.                              *
!             June '91, modified to compute zeta, P, kappa and inte-   *
!             grals for Schwarz inequality in a k2 loop.               *
!             April '92 modified from k2Loop to a separate subroutine  *
!              for estimates of the gradient.                          *
!***********************************************************************

use Real_Spherical, only: ipSph, RSph
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), mZeta, nZeta, nWork2, nHrrMtrx, iShlla, jShllb
real(kind=wp), intent(in) :: CoorM(3,4), Zeta(mZeta), ZInv(mZeta), rKapab(mZeta), P(nZeta,3), rKapcd(mZeta), Q(nZeta,3)
real(kind=wp), intent(out) :: Wrk(nWork2), HMtrx(nHrrMtrx,2)
integer(kind=iwp), intent(out) :: i_Int
integer(kind=iwp) :: iW3, la, lb, mabcd, mabMax, mabMin, mcdMax, mcdMin, ne, nT
real(kind=wp) :: CoorAC(3,2)
logical(kind=iwp), parameter :: NoSpecial = .true.
logical(kind=iwp), external :: EQ
external :: Cff2Dq, ModU2, TERISq, xRys2D
! Statement function to compute canonical index
integer(kind=iwp) :: ixyz, i, nabSz, nElem
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1
nElem(i) = (i+1)*(i+2)/2

la = iAnga(1)
lb = iAnga(2)
#ifdef _DEBUGPRINT_
call RecPrt(' In SchInt: CoorM',' ',CoorM,3,4)
call RecPrt(' In SchInt: P',' ',P,nZeta,3)
call RecPrt(' In SchInt: Q',' ',Q,nZeta,3)
write(u6,*) 'iAnga=',iAnga
#endif

! Compute primitive integrals to be used in the prescreening
! by the Schwarz inequality.

! Compute actual size of [a0|c0] block

mabMin = nabSz(max(la,lb)-1)+1
if (EQ(CoorM(1,1),CoorM(1,2))) mabMin = nabSz(la+lb-1)+1
mabMax = nabSz(la+lb)
mcdMin = nabSz(max(la,lb)-1)+1
if (EQ(CoorM(1,3),CoorM(1,4))) mcdMin = nabSz(la+lb-1)+1
mcdMax = mabMax
mabcd = (mabMax-mabMin+1)*(mcdMax-mcdMin+1)

! Find the proper centers to start of with the angular
! momentum on. If la == lb there will exist an
! ambiguity to which center that angular momentum should
! be accumulated on. In that case we will use A and C of
! the order as defined by the basis functions types.

if (iAnga(1) >= iAnga(2)) then
  call dcopy_(3,CoorM(1,1),1,CoorAC(1,1),1)
  call dcopy_(3,CoorM(1,3),1,CoorAC(1,2),1)
else
  call dcopy_(3,CoorM(1,2),1,CoorAC(1,1),1)
  call dcopy_(3,CoorM(1,4),1,CoorAC(1,2),1)
end if

! Compute [a0|c0], ijkl,a,c

nT = mZeta*1
call Rys(iAnga,nT,Zeta,ZInv,mZeta,Zeta,ZInv,mZeta,P,nZeta,Q,nZeta,rKapab,rKapcd,CoorM,CoorM,CoorAC,mabMin,mabMax,mcdMin,mcdMax, &
         Wrk,nWork2,TERISq,ModU2,Cff2Dq,xRys2D,NoSpecial)

#ifdef _DEBUGPRINT_
call RecPrt(' In SchInt: ijkl,[a0|c0]',' ',Wrk,mZeta,mabcd)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate transformation matrix from intermediate integrals
! to final angular composition, observe that since the 2nd order
! particle matrix is transformed from real spherical harmonics
! to cartesians that the prescreening integrals also are in
! cartesians!

ne = (mabMax-mabMin+1)
call HrrMtrx(HMtrx(:,1),ne,la,lb,CoorM(:,1),CoorM(:,2),.false.,RSph(ipSph(la)),nElem(la),.false.,RSph(ipSph(lb)),nElem(lb))
call HrrMtrx(HMtrx(:,2),ne,la,lb,CoorM(:,3),CoorM(:,4),.false.,RSph(ipSph(la)),nElem(la),.false.,RSph(ipSph(lb)),nElem(lb))
!                                                                      *
!***********************************************************************
!                                                                      *
! Apply a transpose prior to Tnsctl to fake the action of Cntrct.

iW3 = 1+mZeta*mabcd
call DGeTMO(Wrk,mZeta,mZeta,mabcd,Wrk(iW3),mabcd)
call dcopy_(mabcd*mZeta,Wrk(iW3),1,Wrk,1)
call TnsCtl(Wrk,nWork2,mZeta,mabMax,mabMin,mabMax,mabMin,HMtrx(:,1),HMtrx(:,2),la,lb,la,lb,nElem(la),nElem(lb),nElem(la), &
            nElem(lb),iShlla,jShllb,iShlla,jShllb,i_Int)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt(' In SchInt',' ',Wrk(i_Int),mZeta,(nElem(la)*nElem(lb))**2)
#endif

end subroutine SchInt
