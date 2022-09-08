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
! Copyright (C) 1993, Roland Lindh                                     *
!               1993, Per Boussard                                     *
!***********************************************************************

subroutine M1Int( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of the M1 integrals used  *
!         ECP calculations. The operator is the nuclear attraction     *
!         operator times a s-type gaussian function.                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
!             Physics, University of Stockholm, Sweden, October '93.   *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Index_Functions, only: nTri3_Elem1, nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: ia, iAnga(4), ib, iDCRT(0:7), ii, iM1xp, ip, ipAInt, ipIn, ipK, ipPx, ipPy, ipPz, iPrint, ipTmp, ipZ, ipZI, &
                     iRout, iZeta, k, kCnt, kCnttp, kdc, l, lDCRT, LmbdT, mabMax, mabMin, mAInt, mArray, nDCRT, nFlop, nMem, nT
real(kind=wp) :: C(3), Coora(3,4), CoorAC(3,2), Coori(3,4), Fact, Factor, Gmma, PTC2, TC(3), Tmp0, Tmp1
character(len=80) :: Label
logical(kind=iwp) :: NoSpecial
logical(kind=iwp), external :: EQ
external :: Cff2D, Fake, TNAI, XRys2D

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(ZInv)
unused_var(nHer)
unused_var(lOper)
unused_var(iChO)
unused_var(nOrdOp)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 193
iPrint = nPrint(iRout)

if (iPrint >= 49) then
  call RecPrt(' In M1Int: A',' ',A,1,3)
  call RecPrt(' In M1Int: RB',' ',RB,1,3)
  call RecPrt(' In M1Int: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In M1Int: P',' ',P,nZeta,3)
  write(u6,*) ' In M1Int: la,lb=',' ',la,lb
end if

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
mabMin = nTri3_Elem1(max(la,lb)-1)
mabMax = nTri3_Elem1(la+lb)-1
if (EQ(A,RB)) mabMin = nTri3_Elem1(la+lb-1)
mAInt = (mabMax-mabMin+1)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  CoorAC(:,1) = A
else
  CoorAC(:,1) = RB
end if

! Compute FLOP's and size of work array which HRR will use.

call mHrr(la,lb,nFlop,nMem)

! Allocate Scratch for primitives and work area for HRR

ip = 1
ipAInt = ip
k = nTri3_Elem1(la+lb)-nTri3_Elem1(max(la,lb)-1)
ip = ip+nZeta*max(k,nMem)
ipK = ip
ip = ip+nZeta
ipZ = ip
ip = ip+nZeta
ipZI = ip
ip = ip+nZeta
ipPx = ip
ip = ip+nZeta
ipPy = ip
ip = ip+nZeta
ipPz = ip
ip = ip+nZeta
if (ip-1 > nArr*nZeta) then
  call WarningMessage(2,'M1Int: ip-1 > nArr*nZeta')
  write(u6,*) ' nArr,nZeta=',nArr,nZeta
  write(u6,*) ' nMem=',nMem
  call Abend()
end if
ipTmp = ip
mArray = nArr*nZeta-ip+1

Array(ipAInt:ipAInt+nZeta*max(k,nMem)-1) = Zero

! Loop over nuclear centers.

kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
  if (.not. dbsc(kCnttp)%ECP) cycle
  if (dbsc(kCnttp)%nM1 == 0) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)
      Coora(:,1) = A
      Coora(:,2) = RB
      Coori(:,1:2) = Coora(:,1:2)
      if ((.not. EQ(A,RB)) .or. (.not. EQ(A,TC))) then
        Coori(1,1) = Coori(1,1)+One
        !Coora(1,1) = Coora(1,1)+One
      end if
      CoorAC(:,2) = TC
      Coori(:,3) = TC
      Coori(:,4) = TC
      Coora(:,3) = TC
      Coora(:,4) = TC

      do iM1xp=1,dbsc(kCnttp)%nM1
        Gmma = dbsc(kCnttp)%M1xp(iM1xp)

        ! Modify the original basis. Observe that
        ! simplification due to A=B are not valid for the
        ! exponent index, eq. P-A=/=0.

        do iZeta=1,nZeta
          PTC2 = (P(iZeta,1)-TC(1))**2+(P(iZeta,2)-TC(2))**2+(P(iZeta,3)-TC(3))**2
          Tmp0 = Zeta(iZeta)+Gmma
          Tmp1 = exp(-Zeta(iZeta)*Gmma*PTC2/Tmp0)
          Array(ipK+iZeta-1) = rKappa(iZeta)*Tmp1
          Array(ipZ+iZeta-1) = Tmp0
          Array(ipZI+iZeta-1) = One/Tmp0
          Array(ipPx+iZeta-1) = (Zeta(iZeta)*P(iZeta,1)+Gmma*TC(1))/Tmp0
          Array(ipPy+iZeta-1) = (Zeta(iZeta)*P(iZeta,2)+Gmma*TC(2))/Tmp0
          Array(ipPz+iZeta-1) = (Zeta(iZeta)*P(iZeta,3)+Gmma*TC(3))/Tmp0
        end do

        ! Compute integrals with the Rys quadrature.

        nT = nZeta
        NoSpecial = .true.
        call Rys(iAnga,nT,Array(ipZ),Array(ipZI),nZeta,[One],[One],1,Array(ipPx),nZeta,TC,1,Array(ipK),[One],Coori,Coora,CoorAC, &
                 mabmin,mabmax,0,0,Array(ipTmp),mArray,TNAI,Fake,Cff2D,XRys2D,NoSpecial)

        ! Accumulate result for all nuclei. Take the charge on
        ! the center into account.

        Factor = -dbsc(kCnttp)%Charge*dbsc(kCnttp)%M1cf(iM1xp)*Fact
        l = nZeta*mAInt
        Array(ipAInt:ipAInt+l-1) = Array(ipAInt:ipAInt+l-1)+Factor*Array(ipTmp:ipTmp+l-1)
        if (iPrint >= 99) then
          call Recprt(' [a+b,0|A|0] in Array',' ',Array(ipTmp),nZeta,mAInt)
          call RecPrt(' [a+b,0|A|0] in AInt',' ',Array(ipAInt),nZeta,mAInt)
        end if

      end do
    end do
  end do
end do

! Use the HRR to compute the required primitive integrals.

call HRR(la,lb,A,RB,Array(ipAInt),nZeta,nMem,ipIn)
ii = ipAInt+ipIn-1
! Move result
call dcopy_(size(rFinal),Array(ii),1,rFinal,1)

if (iPrint >= 99) then
  write(u6,*) ' Result in M1Int'
  do ia=1,nTri_Elem1(la)
    do ib=1,nTri_Elem1(lb)
      write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
      call RecPrt(Label,' ',rFinal(:,ia,ib,1),nAlpha,nBeta)
    end do
  end do
end if

return

end subroutine M1Int
