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

subroutine PPInt( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: to compute pseudo potential integrals.                       *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Index_util, only: nTri0Elem
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#define _USE_WP_
#include "int_interface.fh"
integer(kind=iwp), parameter :: lproju = 9, imax = 100, kcrs = 1
integer(kind=iwp) :: iA, iAB, iAlpha, iB, iBeta, iCntr, iCnttp, iDCRT(0:7), iExp, intmax, iOff, iOff2, ipA, ipScr, iSh, iStrt, &
                     iZeta, kdc, kSh, kShEnd, kShStr, lcr(kcrs), lDCRT, LmbdT, nArray, ncr(imax), ncrr, nDCRT, &
                     nkcrl(lproju+1,kcrs), nkcru(lproju+1,kcrs), nOp, npot, nPP_S
real(kind=wp) :: C(3), ccr(imax), Fact, TC(3), zcr(imax)
integer(kind=iwp), external :: NrOpr
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
call dcopy_(nZeta*nTri0Elem(la)*nTri0Elem(lb)*nIC,[Zero],0,Final,1)
!                                                                      *
!***********************************************************************
!                                                                      *
nArray = 0
ipScr = 1
intmax = max(nTri0Elem(la),nTri0Elem(lb))
intmax = intmax**2
nArray = nArray+intmax
ipA = ipScr+2*intmax
nArray = nArray+nZeta*nTri0Elem(la)*nTri0Elem(lb)
if (nArray > nZeta*nArr) then
  write(u6,*) 'nArray.gt.nZeta*nArr'
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
    write(u6,*) 'dbsc(iCnttp)%nPP-1.gt.lproju'
    write(u6,*) 'dbsc(iCnttp)%nPP=',nPP_S
    write(u6,*) 'lproju            =',lproju
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
      write(u6,*) ' Pseudo: nPot.gt.imax'
      write(u6,*) '         nPot=',nPot
      write(u6,*) '         imax=',imax
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
  write(u6,*) 'ncr',(ncr(i),i=1,npot)
  write(u6,*) 'zcr',(zcr(i),i=1,npot)
  write(u6,*) 'ccr',(ccr(i),i=1,npot)
  write(u6,*) 'nkcrl',(nkcrl(i,1),i=1,iSh)
  write(u6,*) 'nkcru',(nkcru(i,1),i=1,iSh)
# endif

  do iCntr=1,dbsc(iCnttp)%nCntr
    C(1:3) = dbsc(iCnttp)%Coor(1:3,iCntr)

    ! Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+iCntr)%iStab,dc(kdc+iCntr)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

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

          do iB=1,nTri0Elem(lb)
            do iA=1,nTri0Elem(la)
              iAB = (iB-1)*nTri0Elem(la)+iA
              iOff2 = (iB-1)*nTri0Elem(la)*nZeta+(iA-1)*nZeta+iZeta+ipA-1
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
      call SymAdO(Array(ipA),nZeta,la,lb,nComp,Final,nIC,nOp,lOper,iChO,Fact)
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
  call Unused_real_array(PtChrg)
  call Unused_integer(nGrid)
  call Unused_integer(iAddPot)
end if

end subroutine PPInt
