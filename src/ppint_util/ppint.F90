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

subroutine PPInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,final,nZeta,nIC,nComp,la,lb,A,RB,nHer,Array,nArr,Ccoor,nOrdOp,lOper, &
                 iChO,iStabM,nStabM)
!***********************************************************************
!                                                                      *
! Object: to compute pseudo potential integrals.                       *
!                                                                      *
!***********************************************************************

use Basis_Info
use Center_Info
implicit real*8(A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "oneswi.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC), Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta), rKappa(nZeta), &
       P(nZeta,3), A(3), RB(3), C(3), Array(nZeta*nArr), Ccoor(3), TC(3)
integer lOper(nComp), iStabM(0:nStabM-1), iDCRT(0:7), iChO(nComp)

parameter(lproju=9)
parameter(imax=100,kcrs=1)
real*8 ccr(imax), zcr(imax)
integer nkcrl(lproju+1,kcrs), nkcru(lproju+1,kcrs), lcr(kcrs), ncr(imax)
!                                                                      *
!***********************************************************************
!                                                                      *
! Statement function for Cartesian index

nElem(i) = (i+1)*(i+2)/2
!                                                                      *
!***********************************************************************
!                                                                      *
call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)
!                                                                      *
!***********************************************************************
!                                                                      *
nArray = 0
ipScr = 1
intmax = max(nElem(la),nElem(lb))
intmax = intmax**2
nArray = nArray+intmax
ipA = ipScr+2*intmax
nArray = nArray+nZeta*nElem(la)*nElem(lb)
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

  if (dbsc(iCnttp)%nPP == 0) cycle
  !AOM< Get the "true" (non SO) shells
  nPP_S = 0
  do kSh=dbsc(iCnttp)%iPP,dbsc(iCnttp)%iPP+dbsc(iCnttp)%nPP-1
    ! Skip if a cardholder shell
    if (Shells(kSh)%nExp <= 0) cycle
    ncrr = int(Shells(kSh)%exp(1))
    if (ncrr <= 500) nPP_S = nPP_S+1
  end do
  if (nPP_S == 0) cycle
  !AOM>

  npot = 0
  kShStr = dbsc(iCnttp)%iPP
  kShEnd = kShStr+nPP_S-1
  if (nPP_S-1 > lproju) then
    write(6,*) 'dbsc(iCnttp)%nPP-1.gt.lproju'
    write(6,*) 'dbsc(iCnttp)%nPP=',nPP_S
    write(6,*) 'lproju            =',lproju
    call Abend()
  end if
  lcr(kcrs) = nPP_S-1
  iSh = 0
  iOff = 1
  do kSh=kShStr,kShEnd
    iSh = iSh+1
    nkcrl(iSh,kcrs) = iOff
    nkcru(iSh,kcrs) = iOff+Shells(kSh)%nExp/3-1
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
# ifdef _DEBUGPRINT_
  write(6,*) 'ncr',(ncr(i),i=1,npot)
  write(6,*) 'zcr',(zcr(i),i=1,npot)
  write(6,*) 'ccr',(ccr(i),i=1,npot)
  write(6,*) 'nkcrl',(nkcrl(i,1),i=1,iSh)
  write(6,*) 'nkcru',(nkcru(i,1),i=1,iSh)
# endif

  do iCntr=1,dbsc(iCnttp)%nCntr
    C(1:3) = dbsc(iCnttp)%Coor(1:3,iCntr)

    ! Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+iCntr)%iStab,dc(kdc+iCntr)%nStab,iDCRT,nDCRT)
    Fact = dble(nStabM)/dble(LmbdT)

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)
      !                                                                *
      !*****************************************************************
      !                                                                *
      iZeta = 0
      do iBeta=1,nBeta
        do iAlpha=1,nAlpha
          iZeta = iZeta+1
          call FZero(Array(ipScr),intmax)
          call Pseudo(Alpha(iAlpha),A(1),A(2),A(3),la+1,Beta(iBeta),RB(1),RB(2),RB(3),lb+1,Array(ipScr),intmax,max(la+1,lb+1),ccr, &
                      zcr,nkcrl,nkcru,lcr,ncr,TC(1),TC(2),TC(3),npot)

          do iB=1,nElem(lb)
            do iA=1,nElem(la)
              iAB = (iB-1)*nElem(la)+iA
              iOff2 = (iB-1)*nElem(la)*nZeta+(iA-1)*nZeta+iZeta+ipA-1
              Array(iOff2) = Array(iAB+ipScr-1)
            end do ! iA
          end do   ! iB

        end do     ! iAlpha
      end do       ! iBeta
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Symmetry Adapt

      nOp = NrOpr(iDCRT(lDCRT))
      call SymAdO(Array(ipA),nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,Fact)
    end do ! lDCRT
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do   ! iCntr
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
end if

end subroutine PPInt
