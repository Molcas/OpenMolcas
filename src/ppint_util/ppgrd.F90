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

subroutine PPGrd(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,final,nZeta,la,lb,A,RB,nHer,Array,nArr,Ccoor,nOrdOp,Grad,nGrad,IfGrad,&
                 IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,iStabM,nStabM)

!***********************************************************************
!                                                                      *
! Object: to compute pseudo potential gradient integrals               *
!                                                                      *
!***********************************************************************
use Basis_Info
use Center_Info
use Symmetry_Info, only: iOper
implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "oneswi.fh"
#include "print.fh"
#include "disp.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6), Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta), rKappa(nZeta), &
       P(nZeta,3), A(3), RB(3), C(3), Array(nZeta*nArr), Ccoor(3), TC(3), Grad(nGrad), DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
character*80 Label
integer lOper(nComp), iStabM(0:nStabM-1), iDCRT(0:7), IndGrd(3,2), iuvwx(4), kOp(2), lOp(4), JndGrd(3,4)
logical IfGrad(3,2), JfGrad(3,4)

parameter(lproju=9)
parameter(imax=100,kcrs=1)
real*8 ccr(imax), zcr(imax)
integer nkcrl(lproju+1,kcrs), nkcru(lproju+1,kcrs), lcr(kcrs), ncr(imax)

logical EQ
logical, external :: TF

!                                                                      *
!***********************************************************************
!                                                                      *
! Statement function for Cartesian index
!
nElem(i) = (i+1)*(i+2)/2
!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 122
iPrint = nPrint(iRout)

nDAO = nElem(la)*nElem(lb)
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

nlaplb = max(nElem(la+1),nElem(lb))**2
iplaplb = ipRef+2*nArray
nArray = nArray+nlaplb

! la-1, lb

if (la > 0) then
  nlamlb = max(nElem(la-1),nElem(lb))**2
else
  nlamlb = 0
end if
iplamlb = ipRef+2*nArray
nArray = nArray+nlamlb

! la, lb+1

nlalbp = max(nElem(la),nElem(lb+1))**2
iplalbp = ipRef+2*nArray
nArray = nArray+nlalbp

! la, lb-1

if (lb > 0) then
  nlalbm = max(nElem(la),nElem(lb-1))**2
else
  nlalbm = 0
end if
iplalbm = ipRef+2*nArray
nArray = nArray+nlalbm

if (nArray > nZeta*nArr) then
  write(6,*) 'nArray.gt.nZeta*nArr'
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
kdc = 0
do iCnttp=1,nCnttp
  if (iCnttp > 1) kdc = kdc+dbsc(iCnttp-1)%nCntr
  if (dbsc(iCnttp)%Charge == 0d0) cycle
  if (dbsc(iCnttp)%nPP == 0) cycle
  !AOM< Get the "true" (non SO) shells
  nPP_S = 0
  do kSh=dbsc(iCnttp)%iPP,dbsc(iCnttp)%iPP+dbsc(iCnttp)%nPP-1
    if (Shells(kSh)%nExp <= 0) cycle
    ncrr = int(Shells(kSh)%exp(1))
    if (ncrr <= 500) nPP_S = nPP_S+1
  end do
  if (nPP_S == 0) cycle
  !AOM>

  npot = 0
  kShStr = dbsc(iCnttp)%iPP
  !AOM kShEnd = kShStr+dbsc(iCnttp)%nPP-1
  kShEnd = kShStr+nPP_S-1
  if (dbsc(iCnttp)%nPP-1 > lproju) then
    write(6,*) 'dbsc(iCnttp)%nPP-1.gt.lproju'
    !AOM write(6,*) 'dbsc(iCnttp)%nPP=',dbsc(iCnttp)%nPP
    write(6,*) 'dbsc(iCnttp)%nPP=',nPP_S
    write(6,*) 'lproju          =',lproju
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
      write(6,*) ' Pseudo: nPot.gt.imax'
      write(6,*) '         nPot=',nPot
      write(6,*) '         imax=',imax
      call Abend()
    end if
    iStrt = 1
    do iExp=1,Shells(kSh)%nExp/3
      npot = npot+1
      ncr(npot) = int(Shells(kSh)%exp(iStrt))
      zcr(npot) = Shells(kSh)%exp(iStrt+1)
      ccr(npot) = Shells(kSh)%exp(iStrt+2)
      iStrt = iStrt+3
    end do
  end do
  !write(6,*) 'ncr',(ncr(i),i=1,npot)
  !write(6,*) 'zcr',(zcr(i),i=1,npot)
  !write(6,*) 'ccr',(ccr(i),i=1,npot)
  !write(6,*) 'nkcrl',(nkcrl(i,1),i=1,iSh)
  !write(6,*) 'nkcru',(nkcru(i,1),i=1,iSh)

  do kCnt=1,dbsc(iCnttp)%nCntr
    C(1:3) = dbsc(iCnttp)%Coor(1:3,kCnt)

    ! Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = dble(nStabM)/dble(LmbdT)

    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab
    call ICopy(6,IndGrd,1,JndGrd,1)
    do i=1,3
      do j=1,2
        JfGrad(i,j) = IfGrad(i,j)
      end do
    end do

    nDisp = IndDsp(kdc+kCnt,iIrrep)
    do iCar=0,2
      JfGrad(iCar+1,3) = .false.
      iCmp = 2**iCar
      if (TF(kdc+kCnt,iIrrep,iCmp) .and. .not. dbsc(iCnttp)%pChrg) then
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
          call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+2,Beta(iBeta),RB(1),RB(2),RB(3),lb+1,Array(iplaplb),nlaplb,max(la+2,lb+1),&
                      ccr,zcr,nkcrl,nkcru,lcr,ncr,TC(1),TC(2),TC(3),npot)

          ! la-1, lb

          if (la > 0) then
            call FZero(Array(iplamlb),nlamlb)
            call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la,Beta(iBeta),RB(1),RB(2),RB(3),lb+1,Array(iplamlb),nlamlb,max(la,lb+1),ccr,&
                        zcr,nkcrl,nkcru,lcr,ncr,TC(1),TC(2),TC(3),npot)
          end if

          ! la, lb+1

          call FZero(Array(iplalbp),nlalbp)
          call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+1,Beta(iBeta),RB(1),RB(2),RB(3),lb+2,Array(iplalbp),nlalbp,max(la+1,lb+2),&
                      ccr,zcr,nkcrl,nkcru,lcr,ncr,TC(1),TC(2),TC(3),npot)

          ! la, lb-1

          if (lb > 0) then
            call FZero(Array(iplalbm),nlalbm)
            call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+1,Beta(iBeta),RB(1),RB(2),RB(3),lb,Array(iplalbm),nlalbm,max(la+1,lb),ccr,&
                        zcr,nkcrl,nkcru,lcr,ncr,TC(1),TC(2),TC(3),npot)
          end if

          ! Assemble gradient and store in Final.

          call Assemble_PPGrd(final,nZeta,la,lb,iZeta,Alpha(iAlpha),Beta(iBeta),Array(iplaplb),Array(iplamlb),Array(iplalbp),&
                              Array(iplalbm),JfGrad)

        end do ! iAlpha
      end do   ! iBeta

      !AOM<
      if (abs(Fact-1d0) > 1d-7) call dscal_(nAlpha*nBeta*nElem(la)*nElem(lb)*mGrad,Fact,final,1)
      !AOM>
      if (iPrint >= 99) then
        write(6,*) ' Result in PPGrd'
        write(6,*) JfGrad
        do ia=1,nElem(la)
          do ib=1,nElem(lb)
            do iVec=1,mGrad
              write(Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
              call RecPrt(Label,' ',final(1,ia,ib,iVec),nAlpha,nBeta)
            end do
          end do
        end do
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Distribute contributions to the gradient

      call Distg1X(final,DAO,nZeta,nDAO,mGrad,Grad,nGrad,JfGrad,JndGrd,iuvwx,lOp)

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
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Zeta)
  call Unused_real_array(ZInv)
  call Unused_real_array(rKappa)
  call Unused_real_array(P)
  call Unused_integer(nHer)
  call Unused_real_array(Ccoor)
  call Unused_integer(nOrdOp)
  call Unused_integer_array(lOper)
end if

end subroutine PPGrd
