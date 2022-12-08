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

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Index_Functions, only: nTri_Elem1
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: ia, iab, ib, iDCRT(0:7), iM2xp, ipab, ipAxyz, ipBxyz, ipK, ipPx, ipPy, ipPz, ipQxyz, ipRes, iPrint, ipRxyz, &
                     ipZ, iRout, iZeta, kCnt, kCnttp, kdc, lDCRT, LmbdT, nDCRT, nip
real(kind=wp) :: C(3), Fact, Factor, Gmma, PTC2, TC(3), Tmp0, Tmp1
character(len=80) :: Label
logical(kind=iwp) :: ABeq(3)

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(ZInv)
unused_var(lOper)
unused_var(iChO)
unused_var(PtChrg)
unused_var(iAddPot)

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
nip = nip+nZeta*nComp*nTri_Elem1(la)*nTri_Elem1(lb)
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'M2Int: nip-1 > nArr*nZeta')
  write(u6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(u6,*) ' Abend in M2Int'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In M2Int: A',' ',A,1,3)
  call RecPrt(' In M2Int: RB',' ',RB,1,3)
  call RecPrt(' In M2Int: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In M2Int: Kappa',' ',rKappa,nAlpha,nBeta)
  call RecPrt(' In M2Int: Zeta',' ',Zeta,nAlpha,nBeta)
  call RecPrt(' In M2Int: P',' ',P,nZeta,3)
  write(u6,*) ' In M2Int: la,lb,nHer=',la,lb,nHer
end if

rFinal(:,:,:,:) = Zero

! Loop over nuclear centers

kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
  if (.not. dbsc(kCnttp)%ECP) cycle
  if (dbsc(kCnttp)%nM2 == 0) cycle

  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)

      do iM2xp=1,dbsc(kCnttp)%nM2
        Gmma = dbsc(kCnttp)%M2xp(iM2xp)
        if (iPrint >= 99) write(u6,*) ' Gamma=',Gmma

        ! Modify the original basis.

        do iZeta=1,nZeta
          PTC2 = (P(iZeta,1)-TC(1))**2+(P(iZeta,2)-TC(2))**2+(P(iZeta,3)-TC(3))**2
          Tmp0 = Zeta(iZeta)+Gmma
          Tmp1 = exp(-Zeta(iZeta)*Gmma*PTC2/Tmp0)
          Array(ipK+iZeta-1) = rKappa(iZeta)*Tmp1
          Array(ipZ+iZeta-1) = Tmp0
          Array(ipPx+iZeta-1) = (Zeta(iZeta)*P(iZeta,1)+Gmma*TC(1))/Tmp0
          Array(ipPy+iZeta-1) = (Zeta(iZeta)*P(iZeta,2)+Gmma*TC(2))/Tmp0
          Array(ipPz+iZeta-1) = (Zeta(iZeta)*P(iZeta,3)+Gmma*TC(3))/Tmp0
        end do
        if (iPrint >= 99) then
          write(u6,*) ' The modified basis set'
          call RecPrt(' In M2Int: Kappa',' ',Array(ipK),nAlpha,nBeta)
          call RecPrt(' In M2Int: Zeta',' ',Array(ipZ),nAlpha,nBeta)
          call RecPrt(' In M2Int: P',' ',Array(ipPx),nZeta,3)
        end if

        ! Compute the cartesian values of the basis functions angular part

        ABeq(:) = (A == RB) .and. (A == TC)
        call CrtCmp(Array(ipZ),Array(ipPx),nZeta,A,Array(ipAxyz),la,HerR(iHerR(nHer)),nHer,ABeq)
        call CrtCmp(Array(ipZ),Array(ipPx),nZeta,RB,Array(ipBxyz),lb,HerR(iHerR(nHer)),nHer,ABeq)

        ! Compute the contribution from the multipole moment operator

        ABeq(:) = .false.
        call CrtCmp(Array(ipZ),Array(ipPx),nZeta,TC,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

        ! Compute the cartesian components for the multipole
        ! moment integrals. The integrals are factorized into components.

        call Assmbl(Array(ipQxyz),Array(ipAxyz),la,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb,nZeta,HerW(iHerW(nHer)),nHer)

        ! Combine the cartesian components to the full one electron integral.

        call CmbnMP(Array(ipQxyz),nZeta,la,lb,nOrdOp,Array(ipZ),Array(ipK),Array(ipRes),nComp)
        if (iPrint >= 99) then
          write(u6,*) ' Intermediate result in M2Int'
          do ia=1,nTri_Elem1(la)
            do ib=1,nTri_Elem1(lb)
              iab = (ib-1)*nTri_Elem1(la)+ia
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
        if (iPrint >= 99) write(u6,*) ' Factor=',Factor
        call DaXpY_(size(rFinal),Factor,Array(ipRes),1,rFinal,1)

      end do

    end do
  end do

end do

if (iPrint >= 99) then
  write(u6,*) ' Result in M2Int'
  do ia=1,nTri_Elem1(la)
    do ib=1,nTri_Elem1(lb)
      write(Label,'(A,I2,A,I2,A)') ' rFinal(ia=',ia,',ib=',ib,')'
      call RecPrt(Label,' ',rFinal(:,ia,ib,1),nAlpha,nBeta)
    end do
  end do
end if

return

end subroutine M2Int
