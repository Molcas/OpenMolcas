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

subroutine PrjGrd_mck( &
#                     define _CALLING_
#                     include "grd_mck_interface.fh"
                     )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of ECP integrals.         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
!             Physics, University of Stockholm, Sweden, October 1993.  *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_mck_interface.fh"
integer(kind=iwp) :: iAng, iCnt, iDCRT(0:7), iIrrep, Indx(3,4), ip, ipFA1, ipFA2, ipFB1, ipFB2, ipFin, iShll, iuvwx(4), &
                     JndGrd(3,4,0:7), kCnt, kCnttp, kdc, lDCRT, LmbdT, mOp(4), mvec, nBasisi, nDCRT, nExpi, nt
real(kind=wp) C(3), Dum(1), Fact, TC(3)
logical(kind=iwp) :: DiffCnt, ifg(4), ifhess_dum(3,4,3,4), JfGrad(3,4), tr(4)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ
! Statement function for Cartesian index
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

!                                                                      *
!***********************************************************************
!                                                                      *
iuvwx(1) = iu
iuvwx(2) = iv
mop(1) = nOp(1)
mop(2) = nOp(2)
DiffCnt = IfGrad(iDCar,1) .or. IfGrad(iDCar,2)

#ifdef _DEBUGPRINT_
call RecPrt(' In PrjGrd: A',' ',A,1,3)
call RecPrt(' In PrjGrd: RB',' ',RB,1,3)
call RecPrt(' In PrjGrd: P',' ',P,nZeta,3)
call RecPrt(' In PrjGrd: Alpha',' ',Alpha,nAlpha,1)
call RecPrt(' In PrjGrd: Beta',' ',Beta,nBeta,1)
write(u6,*) ' In PrjGrd: la,lb=',' ',la,lb
write(u6,*) ' In PrjGrd: Diffs=',' ',IfGrad(iDCar,1),IfGrad(iDCar,2)
write(u6,*) ' In PrjGrd: Center=',' ',iDCNT
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
    call LCopy(4,[.false.],0,iFg,1)
    call LCopy(4,[.false.],0,tr,1)
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
      ifg(1) = .true.
      ifg(2) = .true.
      Tr(3) = .true.
      JfGrad(iDCar,1) = .true.
      JfGrad(iDCar,2) = .true.
      do iIrrep=0,nIrrep-1
        jndGrd(iDCar,3,iIrrep) = -IndGrd(iIrrep)
      end do
    end if

    do lDCRT=0,nDCRT-1

      mop(3) = nropr(iDCRT(lDCRT))
      mop(4) = mop(3)
      call OA(iDCRT(lDCRT),C,TC)

      if (EQ(A,RB) .and. EQ(A,TC)) cycle

      do iAng=0,dbsc(kCnttp)%nPrj-1
        iShll = dbsc(kCnttp)%iPrj+iAng
        nExpi = Shells(iShll)%nExp
        nBasisi = Shells(iShll)%nBasis
#       ifdef _DEBUGPRINT_
        write(u6,*) 'nExp(iShll)=',nExpi
        write(u6,*) 'nBasisi=',nBasisi
        write(u6,*) ' iAng=',iAng
        call RecPrt('TC',' ',TC,1,3)
#       endif

        if ((nExpi == 0) .or. (nBasisi == 0)) cycle

        ip = 1

        ipFin = ip
        ip = ip+nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2*6

        ipFA1 = ip
        ip = ip+nAlpha*nExpi*nElem(la)*nElem(iAng)*4

        ipFB1 = ip
        ip = ip+nExpi*nBeta*nElem(iAng)*nElem(lb)*4

        ipFB2 = ip
        ipFA2 = ip
        if (ip >= narr) then
          write(u6,*) 'No mem in prjgrd',ip,narr
          call abend()
        end if

        call dcopy_(nArr,[Zero],0,Array,1)

#       ifdef _DEBUGPRINT_
        call Acore(iang,la,ishll,nordop,TC,A,Array(ip),narr-ip+1,Alpha,nalpha,Array(ipFA1),array(ipFA2),jfgrad(1,1),ifhess_dum,1, &
                   .true.)
#       else
        call Acore(iang,la,ishll,nordop,TC,A,Array(ip),narr-ip+1,Alpha,nalpha,Array(ipFA1),array(ipFA2),jfgrad(1,1),ifhess_dum,1, &
                   .false.)
#       endif
        call LToCore(Array(ipFA1),nalpha,ishll,la,iAng,4)

#       ifdef _DEBUGPRINT_
        call coreB(iang,lb,ishll,nordop,TC,RB,Array(ip),narr-ip+1,Beta,nbeta,Array(ipFB1),array(ipFB2),jfgrad(1,2),ifhess_dum,1, &
                   .true.)
#       else
        call coreB(iang,lb,ishll,nordop,TC,RB,Array(ip),narr-ip+1,Beta,nbeta,Array(ipFB1),array(ipFB2),jfgrad(1,2),ifhess_dum,1, &
                   .false.)
#       endif
        call RToCore(Array(ipFB1),nBeta,ishll,lb,iAng,4)

        call CmbnACB1(Array(ipFA1),Array(ipFB1),Array(ipFin),Fact,nAlpha,nBeta,Dum,nBasisi,la,lb,iang,jfgrad,Dum,.false.,Indx, &
                      mvec,idcar)

        nt = nAlpha*nBeta*nElem(lb)*nElem(la)
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

end subroutine PrjGrd_mck
