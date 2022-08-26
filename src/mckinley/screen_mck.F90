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
! Copyright (C) 1992, Roland Lindh                                     *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Screen_mck(PAO,Scrtch,mPAO,nZeta,nEta,mZeta,mEta,lZeta,lEta,Zeta,ZInv,P,xA,xB,rKA,Data1,IndZ,abmax,Eta,EInv,Q,xG,xD, &
                      rKC,Data2,IndE,cdmax,xpre,iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2,CutInt,PreScr,IndZet,IndEta,ldot)

#define _Old_Code_
#ifdef _Old_Code_

!***********************************************************************
!                                                                      *
! Object: to prescreen the integral derivatives.                       *
!                                                                      *
!   nZeta, nEta : unpartitioned length of primitives.                  *
!                                                                      *
!   mZeta, mEta : section length due to partioning. These are usually  *
!                 equal to nZeta and nEta.                             *
!                                                                      *
!   lZeta, lEta : section length after prescreening.                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             March '92                                                *
!                                                                      *
!             April '92 modified for gradient estimate                 *
!                                                                      *
!             Anders Bernhardsson 1995                                 *
!***********************************************************************

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "ndarray.fh"
integer(kind=iwp) :: mPAO, nZeta, nEta, mZeta, mEta, lZeta, lEta, IndZ(mZeta), IndE(mEta), iphX1, iphY1, iphZ1, iphX2, iphY2, &
                     iphZ2, IndZet(nZeta), IndEta(nEta)
real(kind=wp) :: PAO(mZeta*mEta*mPAO), Scrtch(mZeta*mEta*(1+mPAO*2)), Zeta(nZeta), ZInv(nZeta), P(nZeta,3), xA(nZeta), xB(nZeta), &
                 rKA(nZeta), Data1(nZeta*nDArray), abmax, Eta(nEta), EInv(nEta), Q(nEta,3), xG(nEta), xD(nEta), rKC(nEta), &
                 Data2(nEta*nDArray), cdmax, xpre(mZeta*mEta), CutInt
logical(kind=iwp) :: PreScr, ldot
#ifdef _DEBUGPRINT_
#include "print.fh"
#endif
integer(kind=iwp) :: iEta, ij, ip, ip1, ip2, iPAO, ipOAP, ipPAO, iZE, iZeta, jEta, jPAO, jZeta
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iPrint, iRout
#endif
real(kind=wp) :: abcd, Et, rKAB, rKCD, Zt
integer(kind=iwp), external :: ip_ab, ip_Alpha, ip_Beta, ip_Kappa, ip_PCoor, ip_Z, ip_ZInv

#ifdef _DEBUGPRINT_
iRout = 180
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' Data1',' ',Data1,nZeta,nDArray)
  call RecPrt(' Data2',' ',Data2,nEta,nDArray)
  call RecPrt('2nd order density matrix',' ',PAO,mZeta*mEta,mPAO)
end if
#endif

ip = 1
ipPAO = ip
ip = ip+mZeta*mEta*mPAO

! Compress all indices except zeta

ipOAP = ip
ip = ip+mZeta*mEta*mPAO
if (ldot) call DGetMO(PAO,mZeta,mZeta,mEta*mPAO,Scrtch(ipOAP),mEta*mPAO)

! Prescreen Zeta

lZeta = 0
IndZet(:) = 0
if (PreScr) then
  do iZeta=1,mZeta
    jZeta = IndZ(iZeta)
    IndZet(jZeta) = -lZeta
    abcd = Data1(ip_ab(iZeta,nZeta))*cdMax
    if (abs(abcd) >= CutInt) then
      lZeta = lZeta+1
      IndZet(jZeta) = lZeta
      Zeta(lZeta) = Data1(ip_Z(iZeta,nZeta))
      rKA(lZeta) = Data1(ip_Kappa(iZeta,nZeta))
      P(lZeta,1) = Data1(ip_PCoor(iZeta,nZeta))
      P(lZeta,2) = Data1(ip_PCoor(iZeta+nZeta,nZeta))
      P(lZeta,3) = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))
      xA(lZeta) = Data1(ip_Alpha(iZeta,nZeta,1))
      xB(lZeta) = Data1(ip_Beta(iZeta,nZeta,2))
      ZInv(lZeta) = Data1(ip_ZInv(iZeta,nZeta))
      ip1 = ipOAP+mEta*mPAO*(iZeta-1)
      ip2 = ipPAO+mEta*mPAO*(lZeta-1)
      if (lDot) Scrtch(ip2:ip2+mEta*mPAO-1) = Scrtch(ip1:ip1+mEta*mPAO-1)
    end if
  end do
else
  do iZeta=1,mZeta
    lZeta = lZeta+1
    jZeta = IndZ(iZeta)
    IndZet(jZeta) = lZeta
    Zeta(lZeta) = Data1(ip_Z(iZeta,nZeta))
    rKA(lZeta) = Data1(ip_Kappa(iZeta,nZeta))
    P(lZeta,1) = Data1(ip_PCoor(iZeta,nZeta))
    P(lZeta,2) = Data1(ip_PCoor(iZeta+nZeta,nZeta))
    P(lZeta,3) = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))
    xA(lZeta) = Data1(ip_Alpha(iZeta,nZeta,1))
    xB(lZeta) = Data1(ip_Beta(iZeta,nZeta,2))
    ZInv(lZeta) = Data1(ip_ZInv(iZeta,nZeta))
    ip1 = ipOAP+mEta*mPAO*(iZeta-1)
    ip2 = ipPAO+mEta*mPAO*(lZeta-1)
    if (lDot) Scrtch(ip2:ip2+mEta*mPAO-1) = Scrtch(ip1:ip1+mEta*mPAO-1)
  end do
end if
if (lZeta /= 0) then

  if (iphX1 /= 1) P(1:lZeta,1) = -P(1:lZeta,1)
  if (iphY1 /= 1) P(1:lZeta,2) = -P(1:lZeta,2)
  if (iphZ1 /= 1) P(1:lZeta,3) = -P(1:lZeta,3)

  ! Transpose eta,mPAO,zeta to mPAO,zeta,eta

  if (lDot) call DGetMO(Scrtch(ipPAO),mEta,mEta,mPAO*lZeta,Scrtch(ipOAP),mPAO*lZeta)

  ! Prescreen Eta

  lEta = 0
  IndEta(:) = 0
  if (PreScr) then
    do iEta=1,mEta
      jEta = IndE(iEta)
      IndEta(jEta) = -lEta ! To be removed
      abcd = Data2(ip_ab(iEta,nEta))*abMax
      if (abs(abcd) >= CutInt) then
        lEta = lEta+1
        IndEta(jEta) = lEta
        Eta(lEta) = Data2(ip_Z(iEta,nEta))
        rKC(lEta) = Data2(ip_Kappa(iEta,nEta))
        Q(lEta,1) = Data2(ip_PCoor(iEta,nEta))
        Q(lEta,2) = Data2(ip_PCoor(iEta+nEta,nEta))
        Q(lEta,3) = Data2(ip_PCoor(iEta+2*nEta,nEta))
        xG(lEta) = Data2(ip_Alpha(iEta,nEta,1))
        xD(lEta) = Data2(ip_Beta(iEta,nEta,2))
        EInv(lEta) = Data2(ip_ZInv(iEta,nEta))
        ip1 = ipOAP+mPAO*lZeta*(iEta-1)
        ip2 = ipPAO+mPAO*lZeta*(lEta-1)
        if (ldot) Scrtch(ip2:ip2+lZeta*mPAO-1) = Scrtch(ip1:ip1+lZeta*mPAO-1)
      end if
    end do
  else
    do iEta=1,mEta
      lEta = lEta+1
      jEta = IndE(iEta)
      IndEta(jEta) = lEta
      Eta(lEta) = Data2(ip_Z(iEta,nEta))
      rKC(lEta) = Data2(ip_Kappa(iEta,nEta))
      Q(lEta,1) = Data2(ip_PCoor(iEta,nEta))
      Q(lEta,2) = Data2(ip_PCoor(iEta+nEta,nEta))
      Q(lEta,3) = Data2(ip_PCoor(iEta+2*nEta,nEta))
      xG(lEta) = Data2(ip_Alpha(iEta,nEta,1))
      xD(lEta) = Data2(ip_Beta(iEta,nEta,2))
      EInv(lEta) = Data2(ip_ZInv(iEta,nEta))
      ip1 = ipOAP+mPAO*lZeta*(iEta-1)
      ip2 = ipPAO+mPAO*lZeta*(lEta-1)
      if (ldot) Scrtch(ip2:ip2+lZeta*mPAO-1) = Scrtch(ip1:ip1+lZeta*mPAO-1)
    end do
  end if
  if (lEta /= 0) then

    if (iphX2 /= 1) Q(1:lEta,1) = -Q(1:lEta,1)
    if (iphY2 /= 1) Q(1:lEta,2) = -Q(1:lEta,2)
    if (iphZ2 /= 1) Q(1:lEta,3) = -Q(1:lEta,3)

    ! Transpose mPAO,zeta,eta to zeta,eta,mPAO

    if (ldot) call DGeTMO(Scrtch(ipPAO),mPAO,mPAO,lZeta*lEta,PAO,lZeta*lEta)

  end if
end if

ij = 0
do iEta=1,lEta
  Et = Eta(iEta)
  rKCD = rkC(iEta)
  do iZeta=1,lZeta
    Zt = Zeta(iZeta)
    rKAB = rkA(iZeta)
    ij = ij+1
    xpre(ij) = rKAB*rKCD*sqrt(One/(Zt+Et))
  end do
end do
if (ldot) then
  jPAO = 0
  do iPAO=1,mPAO
    do iZE=0,lZeta*lEta-1
      jPAO = jPAO+1
      PAO(jPAO) = xpre(iZE+1)*PAO(jPAO)
    end do
  end do
end if
#ifdef _DEBUGPRINT_
if (iPrint >= 39) call RecPrt(' PAO',' ',PAO,lZeta*lEta,mPAO)
#endif

return

#else

!***********************************************************************
!                                                                      *
! Object: to prescreen the integral derivatives.                       *
!                                                                      *
!   nZeta, nEta : unpartitioned length of primitives.                  *
!                                                                      *
!   mZeta, mEta : section length due to partioning. These are usually  *
!                 equal to nZeta and nEta.                             *
!                                                                      *
!   lZeta, lEta : section length after prescreening.                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             March '92                                                *
!                                                                      *
!             April '92 modified for gradient estimate                 *
!                                                                      *
!             Anders Bernhardsson 1995                                 *
!***********************************************************************

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "ndarray.fh"
integer(kind=iwp) :: mPAO, nZeta, nEta, mZeta, mEta, lZeta, lEta, IndZ(mZeta), IndE(mEta), iphX1, iphY1, iphZ1, iphX2, iphY2, &
                     iphZ2, IndZet(nZeta), IndEta(nEta)
real(kind=wp) :: PAO(mZeta*mEta*mPAO), Scrtch(mZeta*mEta*(1+mPAO*2)), Zeta(nZeta), ZInv(nZeta), P(nZeta,3), xA(nZeta), xB(nZeta), &
                 rKA(nZeta), Data1(nZeta*nDArray), abmax, Eta(nEta), EInv(nEta), Q(nEta,3), xG(nEta), xD(nEta), rKC(nEta), &
                 Data2(nEta*nDArray), cdmax, xpre(mZeta*mEta), CutInt
logical(kind=iwp) :: PreScr, ldot
#ifdef _DEBUGPRINT_
#include "print.fh"
#endif
integer(kind=iwp) :: iEta, ij, ip, ip1, ip2, iPAO, ipOAP, ipPAO, iZE, iZeta, jEta, jPAO, jZeta
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iPrint, iRout
#endif
real(kind=wp) :: abcd, Et, rKAB, rKCD, Zt
integer(kind=iwp), external :: ip_ab, ip_Alpha, ip_Beta, ip_Kappa, ip_PCoor, ip_Z, ip_ZInv

#ifdef _DEBUGPRINT_
iRout = 180
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' Data1',' ',Data1,nZeta,nDArray)
  call RecPrt(' Data2',' ',Data2,nEta,nDArray)
  call RecPrt('2nd order density matrix',' ',PAO,mZeta*mEta,mPAO)
end if
#endif

ip = 1
ipPAO = ip
ip = ip+mZeta*mEta*mPAO

! Compress all indices except zeta

ipOAP = ip
ip = ip+mZeta*mEta*mPAO
if (ldot) call DGetMO(PAO,mZeta,mZeta,mEta*mPAO,Scrtch(ipOAP),mEta*mPAO)

! Prescreen Zeta

lZeta = 0
if (PreScr) then
  do iZeta=1,mZeta
    jZeta = IndZ(iZeta)
    abcd = Data1(ip_ab(jZeta,nZeta))*cdMax
    if (abs(abcd) >= CutInt) then
      lZeta = lZeta+1
      IndZet(lZeta) = IndZ(iZeta)
      Zeta(lZeta) = Data1(ip_Z(iZeta,nZeta))
      rKA(lZeta) = Data1(ip_Kappa(iZeta,nZeta))
      P(lZeta,1) = Data1(ip_PCoor(iZeta,nZeta))
      P(lZeta,2) = Data1(ip_PCoor(iZeta+nZeta,nZeta))
      P(lZeta,3) = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))
      xA(lZeta) = Data1(ip_Alpha(iZeta,nZeta,1))
      xB(lZeta) = Data1(ip_Beta(iZeta,nZeta,2))
      ZInv(lZeta) = Data1(ip_ZInv(iZeta,nZeta))
      ip1 = ipOAP+mEta*mPAO*(iZeta-1)
      ip2 = ipOAP+mEta*mPAO*(lZeta-1)
      if (lDot) Scrtch(ip2:ip2+mEta*mPAO-1) = Scrtch(ip1:ip1+mEta*mPAO-1)
    end if
  end do
else
  do iZeta=1,mZeta
    lZeta = lZeta+1
    IndZet(lZeta) = IndZ(iZeta)
    Zeta(lZeta) = Data1(ip_Z(iZeta,nZeta))
    rKA(lZeta) = Data1(ip_Kappa(iZeta,nZeta))
    P(lZeta,1) = Data1(ip_PCoor(iZeta,nZeta))
    P(lZeta,2) = Data1(ip_PCoor(iZeta+nZeta,nZeta))
    P(lZeta,3) = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))
    xA(lZeta) = Data1(ip_Alpha(iZeta,nZeta,1))
    xB(lZeta) = Data1(ip_Beta(iZeta,nZeta,2))
    ZInv(lZeta) = Data1(ip_ZInv(iZeta,nZeta))
    ip1 = ipOAP+mEta*mPAO*(iZeta-1)
    ip2 = ipOAP+mEta*mPAO*(lZeta-1)
    if (lDot) Scrtch(ip2:ip2+mEta*mPAO-1) = Scrtch(ip1:ip1+mEta*mPAO-1)
  end do
end if
if (lZeta /= 0) then

  if (iphX1 /= 1) P(1:lZeta,1) = -P(1:lZeta,1)
  if (iphY1 /= 1) P(1:lZeta,2) = -P(1:lZeta,2)
  if (iphZ1 /= 1) P(1:lZeta,3) = -P(1:lZeta,3)

  ! Transpose eta,mPAO,zeta to mPAO,zeta,eta

  if (lDot) call DGetMO(Scrtch(ipOAP),mEta,mEta,mPAO*lZeta,Scrtch(ipPAO),mPAO*lZeta)

  ! Prescreen Eta

  lEta = 0
  if (PreScr) then
    do iEta=1,mEta
      jEta = IndE(iEta)
      abcd = Data2(ip_ab(jEta,nEta))*abMax
      if (abs(abcd) >= CutInt) then
        lEta = lEta+1
        IndEta(lEta) = IndE(iEta)
        Eta(lEta) = Data2(ip_Z(iEta,nEta))
        rKC(lEta) = Data2(ip_Kappa(iEta,nEta))
        Q(lEta,1) = Data2(ip_PCoor(iEta,nEta))
        Q(lEta,2) = Data2(ip_PCoor(iEta+nEta,nEta))
        Q(lEta,3) = Data2(ip_PCoor(iEta+2*nEta,nEta))
        xG(lEta) = Data2(ip_Alpha(iEta,nEta,1))
        xD(lEta) = Data2(ip_Beta(iEta,nEta,2))
        EInv(lEta) = Data2(ip_ZInv(iEta,nEta))
        ip1 = ipPAO+mPAO*lZeta*(iEta-1)
        ip2 = ipPAO+mPAO*lZeta*(lEta-1)
        if (ldot) Scrtch(ip2:ip2+lZeta*mPAO-1) = Scrtch(ip1:ip1+lZeta*mPAO-1)
      end if
    end do
  else
    do iEta=1,mEta
      lEta = lEta+1
      IndEta(lEta) = IndE(iEta)
      Eta(lEta) = Data2(ip_Z(iEta,nEta))
      rKC(lEta) = Data2(ip_Kappa(iEta,nEta))
      Q(lEta,1) = Data2(ip_PCoor(iEta,nEta))
      Q(lEta,2) = Data2(ip_PCoor(iEta+nEta,nEta))
      Q(lEta,3) = Data2(ip_PCoor(iEta+2*nEta,nEta))
      xG(lEta) = Data2(ip_Alpha(iEta,nEta,1))
      xD(lEta) = Data2(ip_Beta(iEta,nEta,2))
      EInv(lEta) = Data2(ip_ZInv(iEta,nEta))
      ip1 = ipPAO+mPAO*lZeta*(iEta-1)
      ip2 = ipPAO+mPAO*lZeta*(lEta-1)
      if (ldot) Scrtch(ip2:ip2+lZeta*mPAO-1) = Scrtch(ip1:ip1+lZeta*mPAO-1)
    end do
  end if
  if (lEta /= 0) then

    if (iphX2 /= 1) Q(1:lEta,1) = -Q(1:lEta,1)
    if (iphY2 /= 1) Q(1:lEta,2) = -Q(1:lEta,2)
    if (iphZ2 /= 1) Q(1:lEta,3) = -Q(1:lEta,3)

    ! Transpose mPAO,zeta,eta to zeta,eta,mPAO

    if (ldot) call DGeTMO(Scrtch(ipPAO),mPAO,mPAO,lZeta*lEta,PAO,lZeta*lEta)

  end if
end if

ij = 0
do iEta=1,lEta
  Et = Eta(iEta)
  rKCD = rkC(iEta)
  do iZeta=1,lZeta
    Zt = Zeta(iZeta)
    rKAB = rkA(iZeta)
    ij = ij+1
    xpre(ij) = rKAB*rKCD*sqrt(One/(Zt+Et))
  end do
end do
if (ldot) then
  jPAO = 0
  do iPAO=1,mPAO
    do iZE=0,lZeta*lEta-1
      jPAO = jPAO+1
      PAO(jPAO) = xpre(iZE+1)*PAO(jPAO)
    end do
  end do
end if
#ifdef _DEBUGPRINT_
if (iPrint >= 39) call RecPrt(' PAO',' ',PAO,lZeta*lEta,mPAO)
#endif

return

#endif

end subroutine Screen_mck
