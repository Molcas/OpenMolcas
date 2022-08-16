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

#define _Old_Code_
#ifdef _Old_Code_

subroutine Screen_mck(PAO,Scrtch,mPAO,nZeta,nEta,mZeta,mEta,lZeta,lEta,Zeta,ZInv,P,xA,xB,rKA,Data1,IndZ,ztmx,abmax,zexpmax,nAlpha, &
                      nBeta,Eta,EInv,Q,xG,xD,rKC,Data2,IndE,etmx,cdmax,eexpmax,nGamma,nDelta,xpre,iphX1,iphY1,iphZ1,iphX2,iphY2, &
                      iphZ2,CutInt,PreScr,IndZet,IndEta,ldot)
!***********************************************************************
!                                                                      *
! Object: to prescreen the integral derivatives.                       *
!                                                                      *
!   nZeta, nEta : unpartioned length of primitives.                    *
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

implicit real*8(A-H,O-Z)
#include "ndarray.fh"
real*8 PAO(mZeta*mEta*mPAO), Scrtch(mZeta*mEta*(1+mPAO*2)), Zeta(nZeta), ZInv(nZeta), P(nZeta,3), Eta(nEta), EInv(nEta), &
       Q(nEta,3), xA(nZeta), xB(nZeta), xG(nEta), xD(nEta), Data1(nZeta*nDArray), Data2(nEta*nDArray), rKA(nZeta), rKC(nEta), &
       xpre(mZeta*mEta)
logical PreScr, ldot
integer IndEta(nEta), IndZet(nZeta), IndZ(mZeta), IndE(mEta)
#include "real.fh"
#ifdef _DEBUGPRINT_
#include "print.fh"
#endif

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
call IZero(IndZet,nZeta)
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
      if (lDot) call dcopy_(mEta*mPAO,Scrtch(ip1),1,Scrtch(ip2),1)
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
    if (lDot) call dcopy_(mEta*mPAO,Scrtch(ip1),1,Scrtch(ip2),1)
  end do
end if
if (lZeta == 0) Go To 999

if (iphX1 /= 1) call DScal_(lZeta,-One,P(1,1),1)
if (iphY1 /= 1) call DScal_(lZeta,-One,P(1,2),1)
if (iphZ1 /= 1) call DScal_(lZeta,-One,P(1,3),1)

! Transpose eta,mPAO,zeta to mPAO,zeta,eta

if (lDot) call DGetMO(Scrtch(ipPAO),mEta,mEta,mPAO*lZeta,Scrtch(ipOAP),mPAO*lZeta)

! Prescreen Eta

lEta = 0
call IZero(IndEta,nEta)
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
      if (ldot) call dcopy_(lZeta*mPAO,Scrtch(ip1),1,Scrtch(ip2),1)
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
    if (ldot) call dcopy_(lZeta*mPAO,Scrtch(ip1),1,Scrtch(ip2),1)
  end do
end if
if (lEta == 0) Go To 999

if (iphX2 /= 1) call DScal_(lEta,-One,Q(1,1),1)
if (iphY2 /= 1) call DScal_(lEta,-One,Q(1,2),1)
if (iphZ2 /= 1) call DScal_(lEta,-One,Q(1,3),1)

! Transpose mPAO,zeta,eta to zeta,eta,mPAO

if (ldot) call DGeTMO(Scrtch(ipPAO),mPAO,mPAO,lZeta*lEta,PAO,lZeta*lEta)

999 continue

ij = 0
do iEta=1,lEta
  Et = Eta(iEta)
  rKCD = rkC(iEta)
  do iZeta=1,lZeta
    Zt = Zeta(iZeta)
    rKAB = rkA(iZeta)
    ij = ij+1
    xpre(ij) = rKAB*rKCD*sqrt(1.0d0/(Zt+Et))
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
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(ztmx)
  call Unused_real(zexpmax)
  call Unused_integer(nAlpha)
  call Unused_integer(nBeta)
  call Unused_real(etmx)
  call Unused_real(eexpmax)
  call Unused_integer(nGamma)
  call Unused_integer(nDelta)
end if

end subroutine Screen_mck

#else

subroutine Screen_mck(PAO,Scrtch,mPAO,nZeta,nEta,mZeta,mEta,lZeta,lEta,Zeta,ZInv,P,xA,xB,rKA,Data1,IndZ,ztmx,abmax,zexpmax,nAlpha, &
                      nBeta,Eta,EInv,Q,xG,xD,rKC,Data2,IndE,etmx,cdmax,eexpmax,nGamma,nDelta,xpre,iphX1,iphY1,iphZ1,iphX2,iphY2, &
                      iphZ2,CutInt,PreScr,IndZet,IndEta,ldot)
!***********************************************************************
!                                                                      *
! Object: to prescreen the integral derivatives.                       *
!                                                                      *
!   nZeta, nEta : unpartioned length of primitives.                    *
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

implicit real*8(A-H,O-Z)
#include "ndarray.fh"
real*8 PAO(mZeta*mEta*mPAO), Scrtch(mZeta*mEta*(1+mPAO*2)), Zeta(nZeta), ZInv(nZeta), P(nZeta,3), Eta(nEta), EInv(nEta), &
       Q(nEta,3), xA(nZeta), xB(nZeta), xG(nEta), xD(nEta), Data1(nZeta*nDArray), Data2(nEta*nDArray), rKA(nZeta), rKC(nEta), &
       xpre(mZeta*mEta)
logical PreScr, ldot
integer IndEta(nEta), IndZet(nZeta), IndZ(mZeta), IndE(mEta)
#include "real.fh"
#ifdef _DEBUGPRINT_
#include "print.fh"
#endif

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
      if (lDot) call dcopy_(mEta*mPAO,Scrtch(ip1),1,Scrtch(ip2),1)
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
    if (lDot) call dcopy_(mEta*mPAO,Scrtch(ip1),1,Scrtch(ip2),1)
  end do
end if
if (lZeta == 0) Go To 999

if (iphX1 /= 1) call DScal_(lZeta,-One,P(1,1),1)
if (iphY1 /= 1) call DScal_(lZeta,-One,P(1,2),1)
if (iphZ1 /= 1) call DScal_(lZeta,-One,P(1,3),1)

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
      if (ldot) call dcopy_(lZeta*mPAO,Scrtch(ip1),1,Scrtch(ip2),1)
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
    if (ldot) call dcopy_(lZeta*mPAO,Scrtch(ip1),1,Scrtch(ip2),1)
  end do
end if
if (lEta == 0) Go To 999

if (iphX2 /= 1) call DScal_(lEta,-One,Q(1,1),1)
if (iphY2 /= 1) call DScal_(lEta,-One,Q(1,2),1)
if (iphZ2 /= 1) call DScal_(lEta,-One,Q(1,3),1)

! Transpose mPAO,zeta,eta to zeta,eta,mPAO

if (ldot) call DGeTMO(Scrtch(ipPAO),mPAO,mPAO,lZeta*lEta,PAO,lZeta*lEta)

999 continue

ij = 0
do iEta=1,lEta
  Et = Eta(iEta)
  rKCD = rkC(iEta)
  do iZeta=1,lZeta
    Zt = Zeta(iZeta)
    rKAB = rkA(iZeta)
    ij = ij+1
    xpre(ij) = rKAB*rKCD*sqrt(1.0d0/(Zt+Et))
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
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(ztmx)
  call Unused_real(zexpmax)
  call Unused_integer(nAlpha)
  call Unused_integer(nBeta)
  call Unused_real(etmx)
  call Unused_real(eexpmax)
  call Unused_integer(nGamma)
  call Unused_integer(nDelta)
end if

end subroutine Screen_mck

#endif
