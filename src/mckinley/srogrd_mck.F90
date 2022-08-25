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

subroutine SroGrd_mck( &
#                     define _CALLING_
#                     include "grd_mck_interface.fh"
                     )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of ECP integrals.         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
!             Physics, University of Stockholm, Sweden, October '93.   *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
#include "grd_mck_interface.fh"
integer(kind=iwp) :: iAng, iCnt, iDCRT(0:7), iIrrep, Indx(3,4), ip, ipFA1, ipFA2, ipFB1, ipFB2, ipFin, ipTmp, iShll, iuvwx(4), &
                     JndGrd(3,4,0:7), kCnt, kCnttp, kdc, lDCRT, LmbdT, mOp(4), mvec, nDCRT, nExpi, nt
real(kind=wp) :: C(3), Fact, TC(3)
logical(kind=iwp) :: DiffCnt, ifg(4), ifhess_dum(3,4,3,4), JfGrad(3,4), tr(4)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ

!                                                                      *
!***********************************************************************
!                                                                      *
iuvwx(1) = iu
iuvwx(2) = iv
mOp(1) = nOp(1)
mOp(2) = nOp(2)

DiffCnt = IfGrad(iDCar,1) .or. IfGrad(iDCar,2)

#ifdef _DEBUGPRINT_
call RecPrt(' In SROGrd: A',' ',A,1,3)
call RecPrt(' In SROGrd: RB',' ',RB,1,3)
call RecPrt(' In SROGrd: P',' ',P,nZeta,3)
call RecPrt(' In SROGrd: Alpha',' ',Alpha,nAlpha,1)
call RecPrt(' In SROGrd: Beta',' ',Beta,nBeta,1)
write(u6,*) ' In SROGrd: la,lb=',' ',la,lb
write(u6,*) ' In SROGrd: Diffs=',' ',IfGrad(iDCar,1),IfGrad(iDCar,2)
write(u6,*) ' In SROGrd: Center=',' ',iDCNT
#endif

kdc = 0
do kCnttp=1,nCnttp

  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
  if (.not. dbsc(kCnttp)%ECP) cycle
  if (dbsc(kCnttp)%nSRO <= 0) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr

    if ((.not. DiffCnt) .and. (kdc+kCnt /= iDCnt)) cycle

    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)
    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab

    call LCopy(12,[.false.],0,JFgrad,1)
    call LCopy(4,[.false.],0,tr,1)
    call LCopy(4,[.false.],0,ifg,1)
    call ICopy(12*nIrrep,[0],0,jndGrd,1)

    do iCnt=1,2
      JfGrad(iDCar,iCnt) = IfGrad(iDCar,iCnt)
    end do

    do ICnt=1,2
      if (ifgrad(idcar,iCnt)) then
        ifg(icnt) = .true.
        do iIrrep=0,nIrrep-1
          jndGrd(iDCar,iCnt,iIrrep) = IndGrd(iIrrep)
        end do
      end if
    end do

    if ((kdc+kCnt) == iDCnt) then
      Tr(3) = .true.
      ifg(1) = .true.
      ifg(2) = .true.
      JfGrad(iDCar,1) = .true.
      JfGrad(iDCar,2) = .true.
      do iIrrep=0,nIrrep-1
        jndGrd(iDCar,3,iIrrep) = -IndGrd(iIrrep)
      end do
    end if

    do lDCRT=0,nDCRT-1

      mOp(3) = nropr(iDCRT(lDCRT))
      mOp(4) = mOp(3)

      call OA(iDCRT(lDCRT),C,TC)

      if (EQ(A,RB) .and. EQ(A,TC)) cycle

      do iAng=0,dbsc(kCnttp)%nSRO-1
        iShll = dbsc(kCnttp)%iSRO+iAng
        nExpi = Shells(iShll)%nExp
#       ifdef _DEBUGPRINT_
        nBasisi = Shells(iShll)%nBasis
        write(u6,*) 'nExpi=',nExpi
        write(u6,*) 'nBasis(iShll)=',nBasisi
        write(u6,*) ' iAng=',iAng
        call RecPrt('TC',' ',TC,1,3)
#       endif

        if (nExpi == 0) cycle

        ip = 1

        ipFin = ip
        ip = ip+nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2*6

        ipTmp = ip
        ip = ip+max(nBeta,nAlpha)*nExpi

        ipFA1 = ip
        ip = ip+nAlpha*nExpi*nTri_Elem1(la)*nTri_Elem1(iAng)*2
        ipFA2 = ip ! Not in use for 1st derivative

        ipFB1 = ip
        ip = ip+nExpi*nBeta*nTri_Elem1(iAng)*nTri_Elem1(lb)*2

        ipFB2 = ip ! Not in use for 1st derivatives

        call dcopy_(nArr,[Zero],0,Array,1)

#ifdef _DEBUGPRINT_
        call Acore(iang,la,ishll,nordop,TC,A,Array(ip),narr-ip+1,Alpha,nalpha,Array(ipFA1),array(ipFA2),jfgrad(1,1),ifhess_dum,1, &
                   .true.)
#else
        call Acore(iang,la,ishll,nordop,TC,A,Array(ip),narr-ip+1,Alpha,nalpha,Array(ipFA1),array(ipFA2),jfgrad(1,1),ifhess_dum,1, &
                   .false.)
#endif
        call LToSph(Array(ipFA1),nalpha,ishll,la,iAng,2)

        call dcopy_(nBeta*nExpi*nTri_Elem1(lb)*nTri_Elem1(iAng)*2,[Zero],0,Array(ipFB1),1)
#ifdef _DEBUGPRINT_
        call coreB(iang,lb,ishll,nordop,TC,RB,Array(ip),narr-ip+1,Beta,nbeta,Array(ipFB1),array(ipFB2),jfgrad(1,2),ifhess_dum,1, &
                   .true.)
#else
        call coreB(iang,lb,ishll,nordop,TC,RB,Array(ip),narr-ip+1,Beta,nbeta,Array(ipFB1),array(ipFB2),jfgrad(1,2),ifhess_dum,1, &
                   .false.)
#endif
        call RToSph(Array(ipFB1),nBeta,ishll,lb,iAng,2)

        call CmbnACB1(Array(ipFA1),Array(ipFB1),Array(ipFin),Fact,nAlpha,nBeta,Shells(iShll)%Akl,nExpi,la,lb,iang,jfgrad, &
                      Array(ipTmp),.true.,Indx,mvec,idcar)

        nt = nAlpha*nBeta*nTri_Elem1(lb)*nTri_Elem1(la)
        call SmAdNa(Array(ipFin),nt,rFinal,mop,loper,JndGrd,iuvwx,JfGrad,Indx,idcar,One,iFG,tr)

      end do
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Zeta)
  call Unused_real_array(ZInv)
  call Unused_real_array(rKappa)
  call Unused_real_array(P)
  call Unused_integer(nHer)
  call Unused_real_array(Ccoor)
  call Unused_logical_array(Trans)
end if

end subroutine SroGrd_mck
