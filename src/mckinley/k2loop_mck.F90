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

subroutine k2Loop_mck(Coor,iAnga,iCmpa,iDCRR,nDCRR,rData,ijCmp,Alpha,nAlpha,Beta,nBeta,Coeff1,iBasn,Coeff2,jBasn,nMemab,Wk002, &
                      m002,Wk003,m003,iStb,jStb)
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

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "ndarray.fh"
integer(kind=iwp), intent(in) :: iAnga(4), iCmpa(4), iDCRR(0:7), nDCRR, ijCmp, nAlpha, nBeta, iBasn, jBasn, nMemab, m002, m003, &
                                 iStb, jStb
real(kind=wp), intent(in) :: Coor(3,2), Alpha(nAlpha), Beta(nBeta), Coeff1(nAlpha,iBasn), Coeff2(nBeta,jBasn)
real(kind=wp), intent(out) :: rData(nAlpha*nBeta*nDArray+nDScalar,nDCRR), Wk002(m002)
real(kind=wp), intent(inout) :: Wk003(m003)
integer(kind=iwp) :: mStb(2), nZeta
real(kind=wp) :: abMax, CoorM(3,4), tmp, Tst, ZtMax
integer(kind=iwp), external :: ip_ab, ip_abMax, ip_Alpha, ip_Beta, ip_EstI, ip_IndZ, ip_Kappa, ip_PCoor, ip_Z, ip_ZetaM, ip_ZInv, &
                               ip_ZtMax
real(kind=wp), external :: EstI

call k2Loop_mck_internal(rData)

! This is to allow type punning without an explicit interface
contains

subroutine k2Loop_mck_internal(rData)

  real(kind=wp), target :: rData(nAlpha*nBeta*nDArray+nDScalar,nDCRR)
  integer(kind=iwp), pointer :: iData(:)
  integer(kind=iwp) :: iZeta, lDCRR

  nZeta = nAlpha*nBeta
  mStb(1) = iStb
  mStb(2) = jStb

  CoorM(1,1) = Coor(1,1)
  CoorM(2,1) = Coor(2,1)
  CoorM(3,1) = Coor(3,1)

  do lDCRR=0,nDCRR-1

    call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
    CoorM(:,3:4) = CoorM(:,1:2)

    ! Compute Zeta, P and kappa.

    call c_f_pointer(c_loc(rData(ip_IndZ(1,nZeta),lDCRR+1)),iData,[nAlpha*nBeta+1])
    call DoZeta(Alpha,nAlpha,Beta,nBeta,CoorM(1,1),CoorM(1,2),rData(ip_PCoor(1,nZeta),lDCRR+1),rData(ip_Z(1,nZeta),lDCRR+1), &
                rData(ip_Kappa(1,nZeta),lDCRR+1),rData(ip_ZInv(1,nZeta),lDCRR+1),rData(ip_Alpha(1,nZeta,1),lDCRR+1), &
                rData(ip_Beta(1,nZeta,2),lDCRR+1),iData)
    nullify(iData)

    call SchInt_mck(CoorM,iAnga,nAlpha,nBeta,nMemab,rData(ip_Z(1,nZeta),lDCRR+1),rData(ip_ZInv(1,nZeta),lDCRR+1), &
                    rData(ip_Kappa(1,nZeta),lDCRR+1),rData(ip_PCoor(1,nZeta),lDCRR+1),nZeta,Wk002,m002,Wk003,m003)

    call PckInt_mck(Wk002,nZeta,ijCmp,rData(ip_ab(1,nZeta),lDCRR+1))
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Estimate the largest contracted integral.

    call c_f_pointer(c_loc(rData(ip_IndZ(1,nZeta),lDCRR+1)),iData,[nAlpha*nBeta+1])
    rData(ip_EstI(nZeta),lDCRR+1) = EstI(rData(ip_Z(1,nZeta),lDCRR+1),rData(ip_Kappa(1,nZeta),lDCRR+1),nAlpha,nBeta,Coeff1,iBasn, &
                                         Coeff2,jBasn,rData(ip_ab(1,nZeta),lDCRR+1),iCmpa(1)*iCmpa(2),Wk002,m002,iData)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Find the largest integral estimate (AO Basis).

    Tst = -One
    do iZeta=0,nZeta-1
      Tst = max(rData(ip_Z(iZeta+1,nZeta),lDCRR+1),Tst)
    end do
    rData(ip_ZetaM(nZeta),lDCRR+1) = tst

    Tst = -One
    ZtMax = Zero
    abMax = Zero
    do iZeta=1,nZeta
      tmp = rData(ip_ab(iZeta,nZeta),lDCRR+1)
      if (Tst < tmp) then
        Tst = tmp
        ZtMax = rData(ip_Z(iZeta,nZeta),lDCRR+1)
        abMax = rData(ip_ab(iZeta,nZeta),lDCRR+1)
      end if
    end do
    rData(ip_ZtMax(nZeta),lDCRR+1) = ZtMax
    rData(ip_abMax(nZeta),lDCRR+1) = abMax
  end do

  return

end subroutine k2Loop_mck_internal

end subroutine k2Loop_mck
