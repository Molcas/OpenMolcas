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
! Copyright (C) 1991,1995,2002, Roland Lindh                           *
!***********************************************************************

subroutine NAInt_GIAO( &
#                     define _CALLING_
#                     include "int_interface.fh"
                     )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of electric field         *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January '91                             *
!                                                                      *
! Modified for explicit code, R. Lindh, February '95.                  *
!                                                                      *
! Modified for GIAOs, R. Lindh, June 2002, Tokyo, Japan.               *
!***********************************************************************

use Basis_Info, only: dbsc, Gaussian_Type, nCnttp, Nuclear_Model, Point_Charge
use Center_Info, only: dc
use Constants, only: Zero, One, Two, Three, Pi, TwoP54
use Definitions, only: wp, iwp

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: iAnga_EF(4), iAnga_NA(4), iComp, iDCRT(0:7), ip3, ipEFInt, ipHRR, ipIn, ipNAInt, iPrint, ipRys, iRout, iZeta, &
                     kab, kCnt, kCnttp, kdc, lab, labcd_EF, labcd_NA, lcd_EF, lcd_NA, lDCRT, llOper, LmbdT, mabMax, mabMin, mArr, &
                     mcdMax_EF, mcdMax_NA, mcdMin_EF, mcdMin_NA, nDCRT, nFLOP, nHRR, nMem, nOp, nT
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), EInv, Eta, Fact, rKappcd, TC(3)
logical(kind=iwp) :: NoSpecial
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ
external :: Fake, MODU2, TERI, TNAI, vCff2D, vRys2D, XCff2D, XRys2D
! Statement function for Cartesian index
integer(kind=iwp) :: nElem, nabSz, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(nHer)
unused_var(Ccoor)
unused_var(PtChrg)
unused_var(iAddPot)
!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 200
iPrint = nPrint(iRout)

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,rFinal,1)

call dcopy_(3,A,1,Coori(1,1),1)
call dcopy_(3,RB,1,Coori(1,2),1)

iAnga_EF(1) = la
iAnga_EF(2) = lb
iAnga_NA(1) = la
iAnga_NA(2) = lb
mabMin = nabSz(max(la,lb)-1)+1
if (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1
mabMax = nabSz(la+lb)
lab = (mabMax-mabMin+1)
kab = nElem(la)*nElem(lb)

iAnga_EF(3) = nOrdOp
iAnga_EF(4) = 0
mcdMin_EF = nabSz(nOrdOp-1)+1
mcdMax_EF = nabSz(nOrdop)
lcd_EF = (mcdMax_EF-mcdMin_EF+1)
labcd_EF = lab*lcd_EF

iAnga_NA(3) = nOrdOp-1
iAnga_NA(4) = 0
mcdMin_NA = nabSz(nOrdOp-2)+1
mcdMax_NA = nabSz(nOrdop-1)
lcd_NA = (mcdMax_NA-mcdMin_NA+1)
labcd_NA = lab*lcd_NA

! Compute Flop's and size of work array which HRR will Use.

call mHRR(la,lb,nFLOP,nMem)
nHRR = max(labcd_EF,labcd_NA,lcd_EF*nMem,lcd_NA*nMem)

! Distribute the work array

mArr = nArr-labcd_EF-nHRR
ipEFInt = 1
ipRys = ipEFInt+nZeta*labcd_EF
ipHRR = ipRys+nZeta*mArr

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if
!
llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do

! Modify Zeta if the two-electron code will be used!

if (Nuclear_Model == Gaussian_Type) then
  do iZeta=1,nZeta
    rKappa(iZeta) = rKappa(iZeta)*(TwoP54/Zeta(iZeta))
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over nuclear centers

kdc = 0
if (nCnttp > 0) kdc = -dbsc(1)%nCntr ! to make sure we start at 0
do kCnttp=1,nCnttp
  kdc = kdc+dbsc(kCnttp)%nCntr
  if (dbsc(kCnttp)%Charge == Zero) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
    if (iPrint >= 99) call RecPrt('C',' ',C,1,3)

    ! Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)
      call dcopy_(3,TC,1,CoorAC(1,2),1)
      call dcopy_(3,TC,1,Coori(1,3),1)
      call dcopy_(3,TC,1,Coori(1,4),1)
      !                                                                *
      !*****************************************************************
      !                                                                *
      !------- Compute integrals with the Rys-Gauss quadrature.        *
      !                                                                *
      !*****************************************************************
      ! 1)                                                             *
      ! Do the EF integrals

      nT = nZeta
      if (Nuclear_Model == Gaussian_Type) then
        NoSpecial = .false.
        Eta = dbsc(kCnttp)%ExpNuc
        EInv = One/Eta
        rKappcd = TwoP54/Eta
        ! Tag on the normalization
        rKappcd = rKappcd*(Eta/Pi)**(Three/Two)
        call Rys(iAnga_EF,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,rKappa,[rKappcd],Coori,Coori,CoorAC,mabMin,mabMax, &
                 mcdMin_EF,mcdMax_EF,Array(ipRys),mArr*nZeta,TERI,MODU2,vCff2D,vRys2D,NoSpecial)
      else if (Nuclear_Model == Point_Charge) then
        NoSpecial = .true.
        call Rys(iAnga_EF,nT,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,rKappa,[One],Coori,Coori,CoorAC,mabMin,mabMax,mcdMin_EF, &
                 mcdMax_EF,Array(ipRys),mArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)
      else
        ! ...more to come...
      end if
      !
      ! The integrals are now ordered as ijkl,e,f
      !
      !  a) Change the order to f,ijkl,e
      !  b) Unfold e to ab, f,ijkl,ab
      !  c) Change the order back to ijkl,ab,f
      !
      ! a)

      call DGetMO(Array(ipRys),nZeta*lab,nZeta*lab,lcd_EF,Array(ipHRR),lcd_EF)

      ! b) Use the HRR to unfold e to ab

      call HRR(la,lb,A,RB,Array(ipHRR),lcd_EF*nZeta,nMem,ipIn)
      ip3 = ipHRR-1+ipIn

      ! c)

      call DGetMO(Array(ip3),lcd_EF,lcd_EF,nZeta*kab,Array(ipEFInt),nZeta*kab)

      ! Stored as nZeta,iElem,jElem,iComp
      !                                                                *
      !*****************************************************************
      ! 2)                                                             *
      ! Do the NA integrals

      if (Nuclear_Model == Gaussian_Type) then
        NoSpecial = .false.
        Eta = dbsc(kCnttp)%ExpNuc
        EInv = One/Eta
        rKappcd = TwoP54/Eta
        ! Tag on the normalization
        rKappcd = rKappcd*(Eta/Pi)**(Three/Two)
        call Rys(iAnga_NA,nT,Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,rKappa,[rKappcd],Coori,Coori,CoorAC,mabMin,mabMax,0,0, &
                 Array(ipRys),mArr*nZeta,TERI,MODU2,vCff2D,vRys2D,NoSpecial)
      else if (Nuclear_Model == Point_Charge) then
        NoSpecial = .true.
        call Rys(iAnga_NA,nT,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,rKappa,[One],Coori,Coori,CoorAC,mabMin,mabMax,0,0, &
                 Array(ipRys),mArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)
      else
        ! ...more to come...
      end if

      ! Use the HRR to compute the required primitive integrals

      call HRR(la,lb,A,RB,Array(ipRys),nZeta,nMem,ipNAInt)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Assemble dV/dB

      call Assemble_dVdB(Array(ipNAInt),Array(ipEFInt),nZeta,la,lb,A,RB,TC)

      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Accumulate contributions

      nOp = NrOpr(iDCRT(lDCRT))
      call SymAdO(Array(ipEFInt),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,-Fact*dbsc(kCnttp)%Charge)

    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (Nuclear_Model == Gaussian_Type) then
  do iZeta=1,nZeta
    rKappa(iZeta) = rKappa(iZeta)*(TwoP54/Zeta(iZeta))
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine NAInt_GIAO
