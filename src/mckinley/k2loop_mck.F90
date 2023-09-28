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
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine k2Loop_mck(Coor,iAnga,iCmpa,iDCRR,nDCRR,k2data, &
                      ijCmp,Alpha,nAlpha,Beta,nBeta,Coeff1,iBasn,Coeff2,jBasn,nMemab,Wk002,m002,Wk003,m003)
!***********************************************************************
!                                                                      *
! Object: to compute zeta, kappa, and P.                               *
!         This is done for all unique pairs of centers                 *
!         generated from the symmetry unique centers A and B.          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             June '91, modified to compute zeta, P, kappa and inte-   *
!             grals for Schwartz inequality in a k2 loop.              *
!             January '92 modified to gradient calculations.           *
!             April '92, modified to use the Cauchy-Schwarz inequality *
!              to estimate the integral derivatives.                   *
!                                                                      *
!             May   '95  modified (simplified) for hessian calculation *
!             By Anders Bernhardsson                                   *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp
use k2_structure, only: k2_type

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), iCmpa(4), iDCRR(0:7), nDCRR, ijCmp, nAlpha, nBeta, iBasn, jBasn, nMemab, m002, m003
real(kind=wp), intent(in) :: Coor(3,2), Alpha(nAlpha), Beta(nBeta), Coeff1(nAlpha,iBasn), Coeff2(nBeta,jBasn)
type(k2_type), intent(inout) :: k2Data(nDCRR)
real(kind=wp), intent(out) :: Wk002(m002)
real(kind=wp), intent(inout) :: Wk003(m003)
integer(kind=iwp) :: iZeta, lDCRR, nZeta
real(kind=wp) :: abMax, CoorM(3,4), tmp, Tst, ZtMax
real(kind=wp), external :: EstI

nZeta = nAlpha*nBeta

CoorM(1,1) = Coor(1,1)
CoorM(2,1) = Coor(2,1)
CoorM(3,1) = Coor(3,1)

do lDCRR=0,nDCRR-1

  call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
  CoorM(:,3:4) = CoorM(:,1:2)

  ! Compute Zeta, P and kappa.

  call DoZeta(Alpha,nAlpha,Beta,nBeta,CoorM(1,1),CoorM(1,2),k2Data(lDCRR+1)%PCoor,k2Data(lDCRR+1)%Zeta,k2Data(lDCRR+1)%Kappa, &
              k2Data(lDCRR+1)%ZInv,k2Data(lDCRR+1)%Alpha,k2Data(lDCRR+1)%Beta,k2Data(lDCRR+1)%IndZ)

  call SchInt_mck(CoorM,iAnga,nAlpha,nBeta,nMemab,k2Data(lDCRR+1)%Zeta,k2Data(lDCRR+1)%ZInv,k2Data(lDCRR+1)%Kappa, &
                  k2Data(lDCRR+1)%PCoor,nZeta,Wk002,m002,Wk003,m003)

  call PckInt_mck(Wk002,nZeta,ijCmp,k2Data(lDCRR+1)%ab)
  !                                                                  *
  !*******************************************************************
  !                                                                  *
  ! Estimate the largest contracted integral.

  k2data(lDCRR+1)%EstI = EstI(k2Data(lDCRR+1)%Zeta,k2Data(lDCRR+1)%Kappa,nAlpha,nBeta,Coeff1,iBasn,Coeff2,jBasn, &
                              k2Data(lDCRR+1)%ab,iCmpa(1)*iCmpa(2),Wk002,m002,k2Data(lDCRR+1)%IndZ)
  !                                                                  *
  !*******************************************************************
  !                                                                  *
  ! Find the largest integral estimate (AO Basis).

  Tst = -One
  do iZeta=1,nZeta
    Tst = max(k2Data(lDCRR+1)%Zeta(iZeta),Tst)
  end do
  k2data(lDCRR+1)%ZetaM = tst

  Tst = -One
  ZtMax = Zero
  abMax = Zero
  do iZeta=1,nZeta
    tmp = k2Data(lDCRR+1)%ab(iZeta)
    if (Tst < tmp) then
      Tst = tmp
      ZtMax = k2Data(lDCRR+1)%Zeta(iZeta)
      abMax = k2Data(lDCRR+1)%ab(iZeta)
    end if
  end do
  k2data(lDCRR+1)%ZtMax = ZtMax
  k2data(lDCRR+1)%abMax = abMax
end do

return

end subroutine k2Loop_mck
