!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine PPGrd( &
#                define _CALLING_
#                include "grd_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: to compute pseudo potential gradient integrals               *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: iOper
use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
integer(kind=iwp), parameter :: lproju = 9, imax = 100, kcrs = 1
integer(kind=iwp) :: i, ia, iAlpha, ib, iBeta, iCar, iCmp, iCnttp, iDCRT(0:7), iExp, iIrrep, iOff, iplalbm, iplalbp, iplamlb, &
                     iplaplb, ipRef, iPrint, iRout, iSh, iStrt, iuvwx(4), iVec, iZeta, JndGrd(3,4), kCnt, kdc, kSh, kShEnd, &
                     kShStr, lcr(kcrs), lDCRT, LmbdT, lOp(4), mGrad, nArray, ncr(imax), ncrr, nDAO, nDCRT, nDisp, &
                     nkcrl(lproju+1,kcrs), nkcru(lproju+1,kcrs), nlalbm, nlalbp, nlamlb, nlaplb, npot, nPP_S
real(kind=wp) :: C(3), ccr(imax), Fact, TC(3), zcr(imax)
character(len=80) :: Label
logical(kind=iwp) :: EQ, JfGrad(3,4)
#include "Molcas.fh"
#include "print.fh"
#include "disp.fh"
logical(kind=iwp), external :: TF

#include "macros.fh"
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(P)
unused_var(nHer)
unused_var(Ccoor(1))
unused_var(nOrdOp)

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 122
iPrint = nPrint(iRout)

nDAO = nTri_Elem1(la)*nTri_Elem1(lb)
iIrrep = 0
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
lOp(1) = iOper(kOp(1))
lOp(2) = iOper(kOp(2))
!                                                                      *
!***********************************************************************
!                                                                      *
! Memory managment

nArray = 0
ipRef = 1

! la+1, lb

nlaplb = max(nTri_Elem1(la+1),nTri_Elem1(lb))**2
iplaplb = ipRef+2*nArray
nArray = nArray+nlaplb

! la-1, lb

if (la > 0) then
  nlamlb = max(nTri_Elem1(la-1),nTri_Elem1(lb))**2
else
  nlamlb = 0
end if
iplamlb = ipRef+2*nArray
nArray = nArray+nlamlb

! la, lb+1

nlalbp = max(nTri_Elem1(la),nTri_Elem1(lb+1))**2
iplalbp = ipRef+2*nArray
nArray = nArray+nlalbp

! la, lb-1

if (lb > 0) then
  nlalbm = max(nTri_Elem1(la),nTri_Elem1(lb-1))**2
else
  nlalbm = 0
end if
iplalbm = ipRef+2*nArray
nArray = nArray+nlalbm

if (nArray > nZeta*nArr) then
  write(u6,*) 'nArray > nZeta*nArr'
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
kdc = 0
do iCnttp=1,nCnttp
  if (iCnttp > 1) kdc = kdc+dbsc(iCnttp-1)%nCntr
  if (dbsc(iCnttp)%Charge == Zero) cycle
  if (dbsc(iCnttp)%nPP == 0) cycle
  !AOM< Get the "true" (non SO) shells
  nPP_S = 0
  do kSh=dbsc(iCnttp)%iPP,dbsc(iCnttp)%iPP+dbsc(iCnttp)%nPP-1
    if (Shells(kSh)%nExp <= 0) cycle
    ncrr = int(Shells(kSh)%Exp(1))
    if (ncrr <= 500) nPP_S = nPP_S+1
  end do
  if (nPP_S == 0) cycle
  !AOM>

  npot = 0
  kShStr = dbsc(iCnttp)%iPP
  !AOM kShEnd = kShStr+dbsc(iCnttp)%nPP-1
  kShEnd = kShStr+nPP_S-1
  if (dbsc(iCnttp)%nPP-1 > lproju) then
    write(u6,*) 'dbsc(iCnttp)%nPP-1 > lproju'
    !AOM write(u6,*) 'dbsc(iCnttp)%nPP   =',dbsc(iCnttp)%nPP
    write(u6,*) 'dbsc(iCnttp)%nPP   =',nPP_S
    write(u6,*) 'lproju             =',lproju
    call Abend()
  end if
  !AOM lcr(kcrs) = dbc(iCnttp)%nPP-1
  lcr(kcrs) = nPP_S-1
  iSh = 0
  iOff = 1
  do kSh=kShStr,kShEnd
    iSh = iSh+1
    nkcrl(iSh,kcrs) = iOff
    nkcru(iSh,kcrs) = iOff+Shells(ksh)%nExp/3-1
    iOff = iOff+Shells(kSh)%nExp/3
    if (nPot > imax) then
      write(u6,*) ' Pseudo: nPot > imax'
      write(u6,*) '         nPot=',nPot
      write(u6,*) '         imax=',imax
      call Abend()
    end if
    iStrt = 1
    do iExp=1,Shells(kSh)%nExp/3
      npot = npot+1
      ncr(npot) = int(Shells(kSh)%Exp(iStrt))
      zcr(npot) = Shells(kSh)%Exp(iStrt+1)
      ccr(npot) = Shells(kSh)%Exp(iStrt+2)
      iStrt = iStrt+3
    end do
  end do
  !write(u6,*) 'ncr',(ncr(i),i=1,npot)
  !write(u6,*) 'zcr',(zcr(i),i=1,npot)
  !write(u6,*) 'ccr',(ccr(i),i=1,npot)
  !write(u6,*) 'nkcrl',(nkcrl(i,1),i=1,iSh)
  !write(u6,*) 'nkcru',(nkcru(i,1),i=1,iSh)

  do kCnt=1,dbsc(iCnttp)%nCntr
    C(1:3) = dbsc(iCnttp)%Coor(1:3,kCnt)

    ! Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab
    JndGrd(:,1:2) = IndGrd
    JfGrad(:,1:2) = IfGrad

    nDisp = IndDsp(kdc+kCnt,iIrrep)
    do iCar=0,2
      JfGrad(iCar+1,3) = .false.
      iCmp = 2**iCar
      if (TF(kdc+kCnt,iIrrep,iCmp) .and. (.not. dbsc(iCnttp)%pChrg)) then
        nDisp = nDisp+1
        if (Direct(nDisp)) then
          JndGrd(iCar+1,1) = abs(JndGrd(iCar+1,1))
          JndGrd(iCar+1,2) = abs(JndGrd(iCar+1,2))
          JndGrd(iCar+1,3) = -nDisp
          JfGrad(iCar+1,1) = .true.
          JfGrad(iCar+1,2) = .true.
        else
          JndGrd(iCar+1,3) = 0
        end if
      else
        JndGrd(iCar+1,3) = 0
      end if
    end do
    JndGrd(:,4) = 0
    JfGrad(:,4) = .false.
    mGrad = 0
    do iCar=1,3
      do i=1,2
        if (JfGrad(iCar,i)) mGrad = mGrad+1
      end do
    end do
    if (mGrad == 0) cycle

    do lDCRT=0,nDCRT-1
      lOp(3) = iDCRT(lDCRT)
      lOp(4) = lOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      if (EQ(A,RB) .and. EQ(A,TC)) cycle
      !                                                                *
      !*****************************************************************
      !                                                                *
      iZeta = 0
      do iBeta=1,nBeta
        do iAlpha=1,nAlpha
          iZeta = iZeta+1

          ! la+1, lb

          call FZero(Array(iplaplb),nlaplb)
          call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+2,Beta(iBeta),RB(1),RB(2),RB(3),lb+1,Array(iplaplb),nlaplb,max(la+2,lb+1), &
                      ccr,zcr,nkcrl,nkcru,lcr,ncr,TC(1),TC(2),TC(3),npot)

          ! la-1, lb

          if (la > 0) then
            call FZero(Array(iplamlb),nlamlb)
            call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la,Beta(iBeta),RB(1),RB(2),RB(3),lb+1,Array(iplamlb),nlamlb,max(la,lb+1),ccr, &
                        zcr,nkcrl,nkcru,lcr,ncr,TC(1),TC(2),TC(3),npot)
          end if

          ! la, lb+1

          call FZero(Array(iplalbp),nlalbp)
          call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+1,Beta(iBeta),RB(1),RB(2),RB(3),lb+2,Array(iplalbp),nlalbp,max(la+1,lb+2), &
                      ccr,zcr,nkcrl,nkcru,lcr,ncr,TC(1),TC(2),TC(3),npot)

          ! la, lb-1

          if (lb > 0) then
            call FZero(Array(iplalbm),nlalbm)
            call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+1,Beta(iBeta),RB(1),RB(2),RB(3),lb,Array(iplalbm),nlalbm,max(la+1,lb),ccr, &
                        zcr,nkcrl,nkcru,lcr,ncr,TC(1),TC(2),TC(3),npot)
          end if

          ! Assemble gradient and store in rFinal.

          call Assemble_PPGrd(rFinal,nZeta,la,lb,iZeta,Alpha(iAlpha),Beta(iBeta),Array(iplaplb),Array(iplamlb),Array(iplalbp), &
                              Array(iplalbm),JfGrad)

        end do ! iAlpha
      end do   ! iBeta

      !AOM<
      if (abs(Fact-One) > 1.0e-7_wp) call dscal_(nAlpha*nBeta*nTri_Elem1(la)*nTri_Elem1(lb)*mGrad,Fact,rFinal,1)
      !AOM>
      if (iPrint >= 99) then
        write(u6,*) ' Result in PPGrd'
        write(u6,*) JfGrad
        do ia=1,nTri_Elem1(la)
          do ib=1,nTri_Elem1(lb)
            do iVec=1,mGrad
              write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
              call RecPrt(Label,' ',rFinal(:,ia,ib,1,iVec),nAlpha,nBeta)
            end do
          end do
        end do
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Distribute contributions to the gradient

      call Distg1X(rFinal,DAO,nZeta,nDAO,mGrad,Grad,nGrad,JfGrad,JndGrd,iuvwx,lOp)

    end do ! lDCRT
!                                                                      *
!***********************************************************************
!                                                                      *
  end do   ! kCnt
end do     ! iCnttp
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine PPGrd
