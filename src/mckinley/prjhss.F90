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

subroutine PrjHss( &
#                 define _CALLING_
#                 include "hss_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of ECP integrals.         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
!             Physics, University of Stockholm, Sweden, October 1993.  *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: iOper
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "hss_interface.fh"
integer(kind=iwp) :: iAng, iDCRT(0:7), ip, ipFA1, ipFA2, ipFB1, ipFB2, ipFin, iShll, iuvwx(4), JndGrd(3,4,0:7), &
                     jndhss(4,3,4,3,0:7), kCnt, kCnttp, kdc, kOp(4), lDCRT, LmbdT, mOp(4), nBasisi, nDCRT, nExpi, nt
real(kind=wp) :: C(3), Coor(3,4), Dum(1), Fact, g2(78), TC(3)
logical(kind=iwp) :: ifg(4), JfGrd(3,4), jfhss(4,3,4,3), tr(4)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ

#include "macros.fh"
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(P)
unused_var(rFinal)
unused_var(nHer)
unused_var(Ccoor)
unused_var(lOper)

iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
mop(1:2) = nOp
kOp(1) = iOper(nOp(1))
kOp(2) = iOper(nOp(2))

Coor(:,1) = A
Coor(:,2) = RB

kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
  if (.not. dbsc(kCnttp)%ECP) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab

    do lDCRT=0,nDCRT-1
      kOp(3) = iDCRT(lDCRT)
      kOp(4) = kOp(3)
      mop(3) = nropr(kop(3))
      mop(4) = mop(3)
      call OA(iDCRT(lDCRT),C,TC)
      Coor(:,3) = TC

      if (EQ(A,RB) .and. EQ(A,TC)) cycle
      call NucInd(coor,kdc+kCnt,ifgrd,ifhss,indgrd,indhss,jfgrd,jfhss,jndgrd,jndhss,tr,ifg)
      do iAng=0,dbsc(kCnttp)%nPrj-1
        iShll = dbsc(kCnttp)%iPrj+iAng
        nExpi = Shells(iShll)%nExp
        nBasisi = Shells(iShll)%nBasis
        if ((nExpi == 0) .or. (nBasisi == 0)) cycle

        ip = 1

        ipFin = ip
        ip = ip+nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*21

        ipFA1 = ip
        ip = ip+nAlpha*nExpi*nTri_Elem1(la)*nTri_Elem1(iAng)*4

        ipFA2 = ip
        ip = ip+nAlpha*nExpi*nTri_Elem1(la)*nTri_Elem1(iAng)*6

        ipFB1 = ip
        ip = ip+nExpi*nBeta*nTri_Elem1(iAng)*nTri_Elem1(lb)*4

        ipFB2 = ip
        ip = ip+nExpi*nBeta*nTri_Elem1(iAng)*nTri_Elem1(lb)*6

        Array(:) = Zero
        ! <a|c>,<a'|c>,<a"|c>
        call Acore(iang,la,ishll,nordop,TC,A,Array(ip),narr-ip+1,Alpha,nalpha,Array(ipFA1),array(ipfa2),jfgrd(1,1),jfhss,2,.false.)
        ! Transform to core orbital
        call LToCore(Array(ipFA1),nalpha,ishll,la,iAng,4)
        call LToCore(Array(ipFA2),nalpha,ishll,la,iAng,6)
        ! <c|b>,<c|b'>,<c|b">
        call coreB(iang,lb,ishll,nordop,TC,RB,Array(ip),narr-ip+1,Beta,nbeta,Array(ipFB1),array(ipfb2),jfgrd(1,2),jfhss,2,.false.)
        ! Transform to core orbital
        call RToCore(Array(ipFB1),nbeta,ishll,lb,iAng,4)
        call RToCore(Array(ipFB2),nbeta,ishll,lb,iAng,6)
        ! Construct complete derivative (contract core)
        call CmbnACB2(Array(ipFa1),Array(ipFa2),Array(ipFb1),Array(ipFb2),Array(ipFin),Fact,nalpha,nbeta,Dum,nBasisi,la,lb,iang, &
                      jfhss,dum,.false.)

        ! contract density
        nt = nZeta*nTri_Elem1(la)*nTri_Elem1(lb)
        g2(:) = Zero
        call dGeMV_('T',nT,21,One,Array(ipFin),nT,DAO,1,Zero,g2,1)

        ! distribute in hessian
        call Distg2(g2,Hess,nHess,JndGrd,JfHss,JndHss,iuvwx,kOp,mop,Tr,IfG)

      end do !iang
    end do !DCR
  end do !cnt
end do !cnttp

return

end subroutine PrjHss
