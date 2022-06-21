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

use Basis_Info
use Center_Info

implicit real*8(A-H,O-Z)
external TNAI, Fake, Cff2D, XRys2D
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables
real*8 C(3), TC(3), CoorAC(3,2), Coori(3,4), Coora(3,4)
character*80 Label
integer iDCRT(0:7), iAnga(4)
logical EQ, NoSpecial
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

iRout = 193
iPrint = nPrint(iRout)

if (iPrint >= 49) then
  call RecPrt(' In M1Int: A',' ',A,1,3)
  call RecPrt(' In M1Int: RB',' ',RB,1,3)
  call RecPrt(' In M1Int: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In M1Int: P',' ',P,nZeta,3)
  write(6,*) ' In M1Int: la,lb=',' ',la,lb
end if

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
mabMin = nabSz(max(la,lb)-1)+1
mabMax = nabSz(la+lb)
if (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1
mAInt = (mabMax-mabMin+1)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if

! Compute FLOP's and size of work array which HRR will use.

call mHrr(la,lb,nFlop,nMem)

! Allocate Scratch for primitives and work area for HRR

ip = 1
ipAInt = ip
k = nabSz(la+lb)-nabSz(max(la,lb)-1)
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
  write(6,*) ' nArr,nZeta=',nArr,nZeta
  write(6,*) ' nMem=',nMem
  call Abend()
end if
ipTmp = ip
mArray = nArr*nZeta-ip+1

call dcopy_(nZeta*max(k,nMem),[Zero],0,Array(ipAInt),1)

! Loop over nuclear centers.

kdc = 0
do kCnttp=1,nCnttp
  if (.not. dbsc(kCnttp)%ECP) Go To 111
  if (dbsc(kCnttp)%nM1 == 0) Go To 111
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = dble(nStabM)/dble(LmbdT)

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)
      call dcopy_(3,A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
      call dcopy_(6,Coora(1,1),1,Coori(1,1),1)
      if ((.not. EQ(A,RB)) .or. (.not. EQ(A,TC))) then
        Coori(1,1) = Coori(1,1)+One
        !Coora(1,1) = Coora(1,1)+One
      end if
      call dcopy_(3,TC,1,CoorAC(1,2),1)
      call dcopy_(3,TC,1,Coori(1,3),1)
      call dcopy_(3,TC,1,Coori(1,4),1)
      call dcopy_(3,TC,1,Coora(1,3),1)
      call dcopy_(3,TC,1,Coora(1,4),1)

      do iM1xp=1,dbsc(kCnttp)%nM1
        Gamma = dbsc(kCnttp)%M1xp(iM1xp)

        ! Modify the original basis. Observe that
        ! simplification due to A=B are not valid for the
        ! exponent index, eq. P-A=/=0.

        do iZeta=1,nZeta
          PTC2 = (P(iZeta,1)-TC(1))**2+(P(iZeta,2)-TC(2))**2+(P(iZeta,3)-TC(3))**2
          Tmp0 = Zeta(iZeta)+Gamma
          Tmp1 = exp(-Zeta(iZeta)*Gamma*PTC2/Tmp0)
          Array(ipK+iZeta-1) = rKappa(iZeta)*Tmp1
          Array(ipZ+iZeta-1) = Tmp0
          Array(ipZI+iZeta-1) = One/Tmp0
          Array(ipPx+iZeta-1) = (Zeta(iZeta)*P(iZeta,1)+Gamma*TC(1))/Tmp0
          Array(ipPy+iZeta-1) = (Zeta(iZeta)*P(iZeta,2)+Gamma*TC(2))/Tmp0
          Array(ipPz+iZeta-1) = (Zeta(iZeta)*P(iZeta,3)+Gamma*TC(3))/Tmp0
        end do

        ! Compute integrals with the Rys quadrature.

        nT = nZeta
        NoSpecial = .true.
        call Rys(iAnga,nT,Array(ipZ),Array(ipZI),nZeta,[One],[One],1,Array(ipPx),nZeta,TC,1,Array(ipK),[One],Coori,Coora,CoorAC, &
                 mabmin,mabmax,0,0,Array(ipTmp),mArray,TNAI,Fake,Cff2D,XRys2D,NoSpecial)

        ! Accumulate result for all nuclei. Take the charge on
        ! the center into account.

        Factor = -dbsc(kCnttp)%Charge*dbsc(kCnttp)%M1cf(iM1xp)*Fact
        call DaXpY_(nZeta*mAInt,Factor,Array(ipTmp),1,Array(ipAInt),1)
        if (iPrint >= 99) then
          call Recprt(' [a+b,0|A|0] in Array',' ',Array(ipTmp),nZeta,mAInt)
          call RecPrt(' [a+b,0|A|0] in AInt',' ',Array(ipAInt),nZeta,mAInt)
        end if

      end do
    end do
  end do
111 kdc = kdc+dbsc(kCnttp)%nCntr
end do

! Use the HRR to compute the required primitive integrals.

call HRR(la,lb,A,RB,Array(ipAInt),nZeta,nMem,ipIn)
ii = ipAInt+ipIn-1
! Move result
call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,Array(ii),1,final,1)

if (iPrint >= 99) then
  write(6,*) ' Result in M1Int'
  do ia=1,nElem(la)
    do ib=1,nElem(lb)
      write(Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
      call RecPrt(Label,' ',final(1,ia,ib,1),nAlpha,nBeta)
    end do
  end do
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Alpha)
  call Unused_real_array(Beta)
  call Unused_real_array(ZInv)
  call Unused_integer(nHer)
  call Unused_integer_array(lOper)
  call Unused_integer_array(iChO)
  call Unused_integer(nOrdOp)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine M1Int
