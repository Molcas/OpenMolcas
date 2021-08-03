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

subroutine NAGrd( &
#                define _CALLING_
#                include "grd_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: to compute the gradient of the nuclear attraction integrals. *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             October '91                                              *
!***********************************************************************

use Basis_Info, only: dbsc, Gaussian_Type, iCnttp_Dummy, nCnttp, Nuclear_Model, Point_Charge
use Center_Info, only: dc
use Constants, only: Zero, One, Two, Three, Pi, TwoP54
use Definitions, only: wp, iwp, u6

implicit none
#define _USE_WP_
#include "grd_interface.fh"
integer(kind=iwp) :: i, iAlpha, iAnga(4), iBeta, iCar, iComp, iDAO, iDCRT(0:7), iIrrep, ipA, ipAOff, ipB, ipBOff, ipDAO, iuvwx(4), &
                     iZeta, j, JndGrd(3,4), jpDAO, kCnt, kCnttp, kdc, lDCRT, LmbdT, lOp(4), mGrad, nArray, nDAO, nDCRT, nDisp, &
                     nip, nRys
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), EInv, Eta, Fact, rKappab, rKappcd, TC(3)
logical(kind=iwp) :: JfGrad(3,4)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: TF
! For normal nuclear attraction
external :: Cff2D, Fake, TNAI1
! Finite nuclei
external :: ModU2, TERI1, vCff2D
#include "Molcas.fh"
#include "disp.fh"
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iPrint, iRout
#include "print.fh"
#endif
! Statement function
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

#ifdef _DEBUGPRINT_
iRout = 150
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  write(u6,*) ' In NAGrd: nArr=',nArr
  nDAO = nElem(la)*nElem(lb)
  call RecPrt('DAO',' ',DAO,nZeta,nDAO)
end if
#endif

nRys = nHer

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
ipDAO = nip
nip = nip+nAlpha*nBeta*nElem(la)*nElem(lb)
if (nip-1 > nZeta*nArr) then
  write(u6,*) ' nip-1 > nZeta*nArr'
  call Abend()
end if
nArray = nZeta*nArr-nip+1

iIrrep = 0
iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
call dcopy_(3,A,1,Coori(1,1),1)
call dcopy_(3,RB,1,Coori(1,2),1)
if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
lOp(1) = kOp(1)
lOp(2) = kOp(2)

ipAOff = ipA
do iBeta=1,nBeta
  call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
  ipAOff = ipAOff+nAlpha
end do

ipBOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
  ipBOff = ipBOff+1
end do

! Modify the density matrix with the prefactor

nDAO = nElem(la)*nElem(lb)
if (Nuclear_Model == Point_Charge) then
  do iDAO=1,nDAO
    do iZeta=1,nZeta
      Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
      DAO(iZeta,iDAO) = Fact*DAO(iZeta,iDAO)
    end do
  end do
end if
!if (iPrint >= 99) call RecPrt('DAO',' ',DAO,nZeta,nDAO)

! Loop over nuclear centers

kdc = 0
do kCnttp=1,nCnttp
  if ((kCnttp /= iCnttp_Dummy) .and. (dbsc(kCnttp)%Charge /= Zero)) then
    do kCnt=1,dbsc(kCnttp)%nCntr
      C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

      call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
      Fact = -dbsc(kCnttp)%Charge*real(nStabM,kind=wp)/real(LmbdT,kind=wp)

      ! Modify the density matrix with prefactors in case of finite nuclei

      if (Nuclear_Model == Gaussian_Type) then
        Eta = dbsc(kCnttp)%ExpNuc
        rKappcd = TwoP54/Eta
        ! Tag on the normalization factor of the nuclear Gaussian
        Fact = Fact*(Eta/Pi)**(Three/Two)
        jpDAO = ipDAO
        do iDAO=1,nDAO
          do iZeta=1,nZeta
            ! On flight modification of Kappa
            rKappab = TwoP54*rKappa(iZeta)/Zeta(iZeta)
            Array(jpDAO) = Fact*DAO(iZeta,iDAO)*rKappab*rKappcd*sqrt(One/(Zeta(iZeta)+Eta))
            jpDAO = jpDAO+1
          end do
        end do
      else if (Nuclear_Model == Point_Charge) then
        call DYaX(nZeta*nDAO,Fact,DAO,1,Array(ipDAO),1)
      else
        write(u6,*) 'NaGrd: Fermi type nuclear distribution not implemented yet!'
        call Abend()
      end if
      iuvwx(3) = dc(kdc+kCnt)%nStab
      iuvwx(4) = dc(kdc+kCnt)%nStab
      call ICopy(6,IndGrd,1,JndGrd,1)
      do i=1,3
        do j=1,2
          JfGrad(i,j) = IfGrad(i,j)
        end do
      end do

      ! Derivatives with respect to the operator is computed via the
      ! translational invariance.

      nDisp = IndDsp(kdc+kCnt,iIrrep)
      do iCar=0,2
        iComp = 2**iCar
        if (TF(kdc+kCnt,iIrrep,iComp) .and. (.not. dbsc(kCnttp)%Frag) .and. (.not. dbsc(kCnttp)%pChrg)) then
          nDisp = nDisp+1
          if (Direct(nDisp)) then
            ! Reset flags for the basis set centers so that we
            ! will explicitly compute the derivatives with
            ! respect to those centers. Activate flag for the
            ! third center so that its derivative will be computed
            ! by the translational invariance.
            JndGrd(iCar+1,1) = abs(JndGrd(iCar+1,1))
            JndGrd(iCar+1,2) = abs(JndGrd(iCar+1,2))
            JndGrd(iCar+1,3) = -nDisp
            JfGrad(iCar+1,1) = .true.
            JfGrad(iCar+1,2) = .true.
            JfGrad(iCar+1,3) = .false.
          else
            JndGrd(iCar+1,3) = 0
            JfGrad(iCar+1,3) = .false.
          end if
        else
          JndGrd(iCar+1,3) = 0
          JfGrad(iCar+1,3) = .false.
        end if
      end do
      ! No derivatives with respect to the fourth center.
      call ICopy(3,[0],0,JndGrd(1,4),1)
      JfGrad(1,4) = .false.
      JfGrad(2,4) = .false.
      JfGrad(3,4) = .false.
      mGrad = 0
      do iCar=1,3
        do i=1,2
          if (JfGrad(iCar,i)) mGrad = mGrad+1
        end do
      end do
      !if (iPrint >= 99) write(u6,*) ' mGrad=',mGrad
      if (mGrad == 0) cycle

      do lDCRT=0,nDCRT-1
        lOp(3) = NrOpr(iDCRT(lDCRT))
        lOp(4) = lOp(3)
        call OA(iDCRT(lDCRT),C,TC)
        call dcopy_(3,TC,1,CoorAC(1,2),1)
        call dcopy_(3,TC,1,Coori(1,3),1)
        call dcopy_(3,TC,1,Coori(1,4),1)

        if (Nuclear_Model == Gaussian_Type) then
          Eta = dbsc(kCnttp)%ExpNuc
          EInv = One/Eta
          call Rysg1(iAnga,nRys,nZeta,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,Coori,Coori, &
                     CoorAC,Array(nip),nArray,TERI1,ModU2,vCff2D,Array(ipDAO),nDAO,Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)
        else if (Nuclear_Model == Point_Charge) then
          call Rysg1(iAnga,nRys,nZeta,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,Coori,Coori, &
                     CoorAC,Array(nip),nArray,TNAI1,Fake,Cff2D,Array(ipDAO),nDAO,Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)
        else
        ! more to come...
        end if

        !call RecPrt('In NaGrd: Grad',' ',Grad,nGrad,1)
      end do
    end do
  end if
  kdc = kdc+dbsc(kCnttp)%nCntr
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Final)
  call Unused_real_array(Ccoor)
  call Unused_integer(nOrdOp)
  call Unused_integer_array(lOper)
end if

end subroutine NAGrd
