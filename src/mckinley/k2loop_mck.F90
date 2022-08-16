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

subroutine k2Loop_mck(Coor,iAnga,iCmpa,iDCRR,nDCRR,data,ijCmp,Alpha,nAlpha,Beta,nBeta,Coeff1,iBasn,Coeff2,jBasn,nMemab,Con,Wk002, &
                      m002,Wk003,m003,Wk004,m004,iStb,jStb)
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

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "ndarray.fh"
#include "real.fh"
#include "disp.fh"
#include "disp2.fh"
real*8 Coor(3,2), CoorM(3,4), Alpha(nAlpha), Beta(nBeta), data(nAlpha*nBeta*nDArray+nDScalar,nDCRR), Coeff1(nAlpha,iBasn), &
       Coeff2(nBeta,jBasn), Wk002(m002), Wk003(m003), Wk004(m004), Con(nAlpha*nBeta)
integer iDCRR(0:7), iAnga(4), iCmpa(4), mStb(2)

call k2Loop_mck_internal(data)

! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Con)
  call Unused_real_array(Wk004)
end if

! This is to allow type punning without an explicit interface
contains

subroutine k2Loop_mck_internal(data)

  use iso_c_binding

  real*8, target :: data(nAlpha*nBeta*nDArray+nDScalar,nDCRR)
  integer, pointer :: iData(:)

  nZeta = nAlpha*nBeta
  mStb(1) = iStb
  mStb(2) = jStb

  CoorM(1,1) = Coor(1,1)
  CoorM(2,1) = Coor(2,1)
  CoorM(3,1) = Coor(3,1)

  do lDCRR=0,nDCRR-1

    call OA(iDCRR(lDCRR),Coor(1:3,2),CoorM(1:3,2))
    call dcopy_(6,CoorM(1,1),1,CoorM(1,3),1)

    ! Compute Zeta, P and kappa.

    call c_f_pointer(c_loc(data(ip_IndZ(1,nZeta),lDCRR+1)),iData,[nAlpha*nBeta+1])
    call DoZeta(Alpha,nAlpha,Beta,nBeta,CoorM(1,1),CoorM(1,2),data(ip_PCoor(1,nZeta),lDCRR+1),data(ip_Z(1,nZeta),lDCRR+1), &
                data(ip_Kappa(1,nZeta),lDCRR+1),data(ip_ZInv(1,nZeta),lDCRR+1),data(ip_Alpha(1,nZeta,1),lDCRR+1), &
                data(ip_Beta(1,nZeta,2),lDCRR+1),iData)
    nullify(iData)

    call SchInt_mck(CoorM,iAnga,iCmpa,nAlpha,nBeta,nMemab,data(ip_Z(1,nZeta),lDCRR+1),data(ip_ZInv(1,nZeta),lDCRR+1), &
                    data(ip_Kappa(1,nZeta),lDCRR+1),data(ip_PCoor(1,nZeta),lDCRR+1),nZeta,Wk002,m002,Wk003,m003)

    call PckInt_mck(Wk002,nZeta,ijCmp,data(ip_ab(1,nZeta),lDCRR+1),data(ip_Z(1,nZeta),lDCRR+1))
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Estimate the largest contracted integral.

    call c_f_pointer(c_loc(data(ip_IndZ(1,nZeta),lDCRR+1)),iData,[nAlpha*nBeta+1])
    data(ip_EstI(nZeta),lDCRR+1) = EstI(data(ip_Z(1,nZeta),lDCRR+1),data(ip_Kappa(1,nZeta),lDCRR+1),nAlpha,nBeta,Coeff1,iBasn, &
                                        Coeff2,jBasn,data(ip_ab(1,nZeta),lDCRR+1),iCmpa(1)*iCmpa(2),Wk002,m002,iData)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Find the largest integral estimate (AO Basis).

    Tst = -One
    do iZeta=0,nZeta-1
      Tst = max(data(ip_Z(iZeta+1,nZeta),lDCRR+1),Tst)
    end do
    data(ip_ZetaM(nZeta),lDCRR+1) = tst

    Tst = -One
    ZtMax = Zero
    abMax = Zero
    do iZeta=1,nZeta
      tmp = data(ip_ab(iZeta,nZeta),lDCRR+1)
      if (Tst < tmp) then
        Tst = tmp
        ZtMax = data(ip_Z(iZeta,nZeta),lDCRR+1)
        abMax = data(ip_ab(iZeta,nZeta),lDCRR+1)
      end if
    end do
    data(ip_ZtMax(nZeta),lDCRR+1) = ZtMax
    data(ip_abMax(nZeta),lDCRR+1) = abMax
  end do

  return

end subroutine k2Loop_mck_internal

end subroutine k2Loop_mck
