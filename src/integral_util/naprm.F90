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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NAPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,final,nZeta,nComp,la,lb,A,RB,nRys,Array,nArr,CCoor,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January 1991                            *
!***********************************************************************

use Basis_Info, only: Nuclear_Model, Gaussian_Type, mGaussian_Type, DBSC, Point_Charge
use Constants, only: Zero, One, Two, Three, Pi, TwoP54

implicit none
! Used for normal nuclear attraction integrals
external TNAI, Fake, XCff2D, XRys2D
! Used for finite nuclei
external TERI, ModU2, vCff2D, vRys2D
#include "oneswi.fh"
integer nZeta, la, lb, nComp, nAlpha, nBeta, nArr, nRys, nOrdOp
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp), Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta), rKappa(nZeta), &
       P(nZeta,3), A(3), RB(3), CCoor(3,2), Array(nZeta*nArr)
real*8 C(3), Coora(3,4), Coori(3,4), CoorAC(3,2)
logical EQ, NoSpecial
integer iAnga(4)
integer ixyz, nElem, nabSz, lc, ld, mabMin, mabMax, iZeta, iCnttp, nT, mcdMin, mcdMax, ipOff, mArr, ipIn, nFlop, nMem
real*8 Q_Nuc, Eta, EInv, rKappCD
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

call FZero(final,nZeta*nElem(la)*nElem(lb)*nComp)

lc = 0
ld = 0
iAnga(1) = la
iAnga(2) = lb
iAnga(3) = lc
iAnga(4) = ld
call dcopy_(3,A,1,Coora(1,1),1)
call dcopy_(3,RB,1,Coora(1,2),1)
call dcopy_(2*3,Coora,1,Coori,1)
mabMin = nabSz(max(la,lb)-1)+1
mabMax = nabSz(la+lb)
if (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1

! Compute FLOPs and size of work array which Hrr will use.

call mHrr(la,lb,nFLOP,nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if

! Modify Zeta if the two-electron code will be used!

if ((Nuclear_Model == Gaussian_Type) .or. (Nuclear_Model == mGaussian_Type)) then
  do iZeta=1,nZeta
    rKappa(iZeta) = rKappa(iZeta)*(TwoP54/Zeta(iZeta))
  end do
end if

iCnttp = int(CCoor(1,2))
Q_Nuc = dbsc(iCnttp)%Charge

if (Q_Nuc == Zero) Go To 111
call dcopy_(3,CCoor,1,C,1)
#ifdef _DEBUGPRINT_
call RecPrt('C',' ',C,1,3)
#endif

call DCopy_(3,C,1,CoorAC(1,2),1)
call DCopy_(3,C,1,Coori(1,3),1)
call DCopy_(3,C,1,Coori(1,4),1)
call DCopy_(3,C,1,Coora(1,3),1)
call DCopy_(3,C,1,Coora(1,4),1)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute integrals with the Rys quadrature.
!                                                                      *
!***********************************************************************
!                                                                      *
nT = nZeta
!                                                                      *
!***********************************************************************
!                                                                      *
if (Nuclear_Model == Gaussian_Type) then

  ! Gaussian nuclear charge distribution

  NoSpecial = .false.
  Eta = dbsc(iCnttp)%ExpNuc
  EInv = One/Eta
  rKappcd = TwoP54/Eta
  ! Tag on the normalization
  rKappcd = rKappcd*(Eta/Pi)**(Three/Two)
  ! s-type function
  mcdMin = 0
  mcdMax = 0
  call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,C,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabmin,mabmax,mcdMin,mcdMax, &
           Array,nArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (Nuclear_Model == mGaussian_Type) then

  ! Modified Gaussian nuclear charge distribution

  NoSpecial = .false.
  Eta = dbsc(iCnttp)%ExpNuc
  EInv = One/Eta
  rKappcd = TwoP54/Eta
  ! Tag on the normalization
  rKappcd = rKappcd*(Eta/Pi)**(Three/Two)/(One+Three*dbsc(iCnttp)%w_mGauss/(Two*Eta))
  ! s type function
  mcdMin = 0
  mcdMax = 0
  call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,C,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabmin,mabmax,mcdMin,mcdMax, &
           Array,nArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)

  ! d type function w*(x**2+y**2+z**2)
  if (dbsc(iCnttp)%w_mGauss > Zero) then
    rKappcd = rKappcd*dbsc(iCnttp)%w_mGauss
    iAnga(3) = 2
    mcdMin = nabSz(2+ld-1)+1
    mcdMax = nabSz(2+ld)
    ! tweak the pointers
    ipOff = 1+nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
    mArr = nArr-(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
    call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,C,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabMin,mabMax,mcdMin,mcdMax, &
             Array(ipOff),mArr*nZeta,TERI,ModU2,vCff2D,vRys2D,NoSpecial)
    iAnga(3) = 0

    ! Add the s and d contributions together!

    call Assemble_mGauss(Array,Array(ipOff),nZeta*(mabMax-mabMin+1))
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (Nuclear_Model == Point_Charge) then

  ! Point-like nuclear charge distribution

  NoSpecial = .true.
  Eta = One
  EInv = One
  rKappcd = One
  mcdMin = 0
  mcdMax = 0
  call Rys(iAnga,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,C,1,rKappa,[rKappcd],Coori,Coora,CoorAC,mabMin,mabMax,mcdMin,mcdMax, &
           Array,nArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Use the HRR to compute the required primitive integrals.

call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)

call DCopy_(nZeta*nElem(la)*nElem(lb)*nComp,Array(ipIn),1,final,1)
call DScal_(nZeta*nElem(la)*nElem(lb)*nComp,-Q_Nuc,final,1)

111 continue

if ((Nuclear_Model == Gaussian_Type) .or. (Nuclear_Model == mGaussian_Type)) then
  do iZeta=1,nZeta
    rKappa(iZeta) = rKappa(iZeta)/(TwoP54/Zeta(iZeta))
  end do
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Alpha)
  call Unused_real_array(Beta)
  call Unused_integer(nRys)
  call Unused_integer(nOrdOp)
end if

end subroutine NAPrm
