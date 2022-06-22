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

subroutine M2Int( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of M2 integrals used in   *
!         ECP calculations. The operator is a s-type gaussian          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
!             Physics, University of Stockholm, Sweden, October '93.   *
!***********************************************************************

use Basis_Info
use Center_Info
use Her_RW, only: HerR, HerW, iHerR, iHerW

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables.
real*8 TC(3), C(3)
character*80 Label
logical ABeq(3)
integer iDCRT(0:7)
! Statement function for Cartesian index
nElem(k) = (k+1)*(k+2)/2

iRout = 122
iPrint = nPrint(iRout)

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+1)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+1)
ipRxyz = nip
nip = nip+nZeta*3*nHer
ipQxyz = nip
nip = nip+nZeta*3*(la+1)*(lb+1)
ipK = nip
nip = nip+nZeta
ipZ = nip
nip = nip+nZeta
ipPx = nip
nip = nip+nZeta
ipPy = nip
nip = nip+nZeta
ipPz = nip
nip = nip+nZeta
ipRes = nip
nip = nip+nZeta*nComp*nElem(la)*nElem(lb)
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'M2Int: nip-1 > nArr*nZeta')
  write(6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(6,*) ' Abend in M2Int'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In M2Int: A',' ',A,1,3)
  call RecPrt(' In M2Int: RB',' ',RB,1,3)
  call RecPrt(' In M2Int: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In M2Int: Kappa',' ',rKappa,nAlpha,nBeta)
  call RecPrt(' In M2Int: Zeta',' ',Zeta,nAlpha,nBeta)
  call RecPrt(' In M2Int: P',' ',P,nZeta,3)
  write(6,*) ' In M2Int: la,lb,nHer=',la,lb,nHer
end if

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

! Loop over nuclear centers

kdc = 0
if (nCnttp > 0) kdc = -dbsc(1)%nCntr ! to make sure we start at 0
do kCnttp=1,nCnttp
  kdc = kdc+dbsc(kCnttp)%nCntr
  if (.not. dbsc(kCnttp)%ECP) cycle
  if (dbsc(kCnttp)%nM2 == 0) cycle

  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = dble(nStabM)/dble(LmbdT)

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)

      do iM2xp=1,dbsc(kCnttp)%nM2
        Gamma = dbsc(kCnttp)%M2xp(iM2xp)
        if (iPrint >= 99) write(6,*) ' Gamma=',Gamma

        ! Modify the original basis.

        do iZeta=1,nZeta
          PTC2 = (P(iZeta,1)-TC(1))**2+(P(iZeta,2)-TC(2))**2+(P(iZeta,3)-TC(3))**2
          Tmp0 = Zeta(iZeta)+Gamma
          Tmp1 = exp(-Zeta(iZeta)*Gamma*PTC2/Tmp0)
          Array(ipK+iZeta-1) = rKappa(iZeta)*Tmp1
          Array(ipZ+iZeta-1) = Tmp0
          Array(ipPx+iZeta-1) = (Zeta(iZeta)*P(iZeta,1)+Gamma*TC(1))/Tmp0
          Array(ipPy+iZeta-1) = (Zeta(iZeta)*P(iZeta,2)+Gamma*TC(2))/Tmp0
          Array(ipPz+iZeta-1) = (Zeta(iZeta)*P(iZeta,3)+Gamma*TC(3))/Tmp0
        end do
        if (iPrint >= 99) then
          write(6,*) ' The modified basis set'
          call RecPrt(' In M2Int: Kappa',' ',Array(ipK),nAlpha,nBeta)
          call RecPrt(' In M2Int: Zeta',' ',Array(ipZ),nAlpha,nBeta)
          call RecPrt(' In M2Int: P',' ',Array(ipPx),nZeta,3)
        end if

        ! Compute the cartesian values of the basis functions angular part

        ABeq(1) = (A(1) == RB(1)) .and. (A(1) == TC(1))
        ABeq(2) = (A(2) == RB(2)) .and. (A(2) == TC(2))
        ABeq(3) = (A(3) == RB(3)) .and. (A(3) == TC(3))
        call CrtCmp(Array(ipZ),Array(ipPx),nZeta,A,Array(ipAxyz),la,HerR(iHerR(nHer)),nHer,ABeq)
        call CrtCmp(Array(ipZ),Array(ipPx),nZeta,RB,Array(ipBxyz),lb,HerR(iHerR(nHer)),nHer,ABeq)

        ! Compute the contribution from the multipole moment operator

        ABeq(1) = .false.
        ABeq(2) = .false.
        ABeq(3) = .false.
        call CrtCmp(Array(ipZ),Array(ipPx),nZeta,TC,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

        ! Compute the cartesian components for the multipole
        ! moment integrals. The integrals are factorized into components.

        call Assmbl(Array(ipQxyz),Array(ipAxyz),la,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb,nZeta,HerW(iHerW(nHer)),nHer)

        ! Combine the cartesian components to the full one electron integral.

        call CmbnMP(Array(ipQxyz),nZeta,la,lb,nOrdOp,Array(ipZ),Array(ipK),Array(ipRes),nComp)
        if (iPrint >= 99) then
          write(6,*) ' Intermediate result in M2Int'
          do ia=1,nElem(la)
            do ib=1,nElem(lb)
              iab = (ib-1)*nElem(la)+ia
              ipab = (iab-1)*nZeta+ipRes
              write(Label,'(A,I2,A,I2,A)') ' Array(',ia,',',ib,')'
              if (nComp /= 1) then
                call RecPrt(Label,' ',Array(ipab),nZeta,nComp)
              else
                call RecPrt(Label,' ',Array(ipab),nAlpha,nBeta)
              end if
            end do
          end do
        end if

        ! Multiply result by Zeff*Const

        Factor = -dbsc(kCnttp)%Charge*dbsc(kCnttp)%M2cf(iM2xp)*Fact
        if (iPrint >= 99) write(6,*) ' Factor=',Factor
        call DaXpY_(nZeta*nElem(la)*nElem(lb)*nIC,Factor,Array(ipRes),1,final,1)

      end do

    end do
  end do

end do

if (iPrint >= 99) then
  write(6,*) ' Result in M2Int'
  do ia=1,nElem(la)
    do ib=1,nElem(lb)
      write(Label,'(A,I2,A,I2,A)') ' Final(ia=',ia,',ib=',ib,')'
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
  call Unused_integer_array(lOper)
  call Unused_integer_array(iChO)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine M2Int
