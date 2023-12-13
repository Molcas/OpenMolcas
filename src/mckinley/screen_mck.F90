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

!#define _DEBUGPRINT_
subroutine Screen_mck(iOffZ,iOffE,PAO,Scrtch,mPAO,nZeta,nEta,mZeta,mEta,lZeta,lEta, &
                      k2Data1,k2Data2, &
                      Zeta,ZInv,P,xA,xB,rKA,Eta,EInv,Q,xG,xD,rKC, &
                      xpre,iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2,CutInt,PreScr,IndZet,IndEta,ldot)

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

use k2_structure, only: k2_type
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iOffZ, iOffE, mPAO, nZeta, nEta, mZeta, mEta, iphX1, iphY1, iphZ1, iphX2, iphY2, iphZ2
real(kind=wp), intent(inout) :: PAO(mZeta*mEta*mPAO)
real(kind=wp), intent(out) :: Scrtch(mZeta*mEta*(1+mPAO*2)), Zeta(nZeta), ZInv(nZeta), P(nZeta,3), xA(nZeta), xB(nZeta), &
                              rKA(nZeta), Eta(nEta), EInv(nEta), Q(nEta,3), xG(nEta), xD(nEta), rKC(nEta), xpre(mZeta*mEta)
integer(kind=iwp), intent(out) :: lZeta, lEta, IndZet(nZeta), IndEta(nEta)
type(k2_type), intent(in) :: k2Data1, k2Data2
real(kind=wp), intent(in) :: CutInt
logical(kind=iwp), intent(in) :: PreScr, ldot
integer(kind=iwp) :: iEta, ij, ip, ip1, ip2, iPAO, ipOAP, ipPAO, iZeta, jEta, jPAO, jZeta
real(kind=wp) :: abcd, abMax, cdMax

#ifdef _DEBUGPRINT_
call RecPrt('2nd order density matrix',' ',PAO,mZeta*mEta,mPAO)
#endif

abMax = k2Data1%abMax
cdMax = k2Data2%abMax

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
    jZeta = k2Data1%IndZ(iOffZ+iZeta)
    IndZet(jZeta) = -lZeta
    abcd = k2Data1%ab(iOffZ+iZeta)*cdMax
    if (abs(abcd) >= CutInt) then
      lZeta = lZeta+1
      IndZet(jZeta) = lZeta
      Zeta(lZeta) = k2Data1%Zeta(iOffZ+iZeta)
      rKA(lZeta) = k2Data1%Kappa(iOffZ+iZeta)
      P(lZeta,1) = k2Data1%PCoor(iOffZ+iZeta,1)
      P(lZeta,2) = k2Data1%PCoor(iOffZ+iZeta,2)
      P(lZeta,3) = k2Data1%PCoor(iOffZ+iZeta,3)
      xA(lZeta) = k2Data1%Alpha(iOffZ+iZeta)
      xB(lZeta) = k2Data1%Beta(iOffZ+iZeta)
      ZInv(lZeta) = k2Data1%ZInv(iOffZ+iZeta)
      ip1 = ipOAP+mEta*mPAO*(iZeta-1)
      ip2 = ipPAO+mEta*mPAO*(lZeta-1)
      if (ldot) Scrtch(ip2:ip2+mEta*mPAO-1) = Scrtch(ip1:ip1+mEta*mPAO-1)
    end if
  end do
else
  do iZeta=1,mZeta
    lZeta = lZeta+1
    jZeta = k2Data1%IndZ(iOffZ+iZeta)
    IndZet(jZeta) = lZeta
    Zeta(lZeta) = k2Data1%Zeta(iOffZ+iZeta)
    rKA(lZeta) = k2Data1%Kappa(iOffZ+iZeta)
    P(lZeta,1) = k2Data1%PCoor(iOffZ+iZeta,1)
    P(lZeta,2) = k2Data1%PCoor(iOffZ+iZeta,2)
    P(lZeta,3) = k2Data1%PCoor(iOffZ+iZeta,3)
    xA(lZeta) = k2Data1%Alpha(iOffZ+iZeta)
    xB(lZeta) = k2Data1%Beta(iOffZ+iZeta)
    ZInv(lZeta) = k2Data1%ZInv(iOffZ+iZeta)
    ip1 = ipOAP+mEta*mPAO*(iZeta-1)
    ip2 = ipPAO+mEta*mPAO*(lZeta-1)
    if (ldot) Scrtch(ip2:ip2+mEta*mPAO-1) = Scrtch(ip1:ip1+mEta*mPAO-1)
  end do
end if
if (lZeta /= 0) then

  if (iphX1 /= 1) P(1:lZeta,1) = -P(1:lZeta,1)
  if (iphY1 /= 1) P(1:lZeta,2) = -P(1:lZeta,2)
  if (iphZ1 /= 1) P(1:lZeta,3) = -P(1:lZeta,3)

  ! Transpose eta,mPAO,zeta to mPAO,zeta,eta

  if (ldot) call DGetMO(Scrtch(ipPAO),mEta,mEta,mPAO*lZeta,Scrtch(ipOAP),mPAO*lZeta)

  ! Prescreen Eta

  lEta = 0
  IndEta(:) = 0
  if (PreScr) then
    do iEta=1,mEta
      jEta = k2Data2%IndZ(iOffE+iEta)
      IndEta(jEta) = -lEta ! To be removed
      abcd = k2Data2%ab(iOffE+iEta)*abMax
      if (abs(abcd) >= CutInt) then
        lEta = lEta+1
        IndEta(jEta) = lEta
        Eta(lEta) = k2Data2%Zeta(iOffE+iEta)
        rKC(lEta) = k2Data2%Kappa(iOffE+iEta)
        Q(lEta,1) = k2Data2%PCoor(iOffE+iEta,1)
        Q(lEta,2) = k2Data2%PCoor(iOffE+iEta,2)
        Q(lEta,3) = k2Data2%PCoor(iOffE+iEta,3)
        xG(lEta) = k2Data2%Alpha(iOffE+iEta)
        xD(lEta) = k2Data2%Beta(iOffE+iEta)
        EInv(lEta) = k2Data2%ZInv(iOffE+iEta)
        ip1 = ipOAP+mPAO*lZeta*(iEta-1)
        ip2 = ipPAO+mPAO*lZeta*(lEta-1)
        if (ldot) Scrtch(ip2:ip2+lZeta*mPAO-1) = Scrtch(ip1:ip1+lZeta*mPAO-1)
      end if
    end do
  else
    do iEta=1,mEta
      lEta = lEta+1
      jEta = k2Data2%IndZ(iOffE+iEta)
      IndEta(jEta) = lEta
      Eta(lEta) = k2Data2%Zeta(iOffE+iEta)
      rKC(lEta) = k2Data2%Kappa(iOffE+iEta)
      Q(lEta,1) = k2Data2%PCoor(iOffE+iEta,1)
      Q(lEta,2) = k2Data2%PCoor(iOffE+iEta,2)
      Q(lEta,3) = k2Data2%PCoor(iOffE+iEta,3)
      xG(lEta) = k2Data2%Alpha(iOffE+iEta)
      xD(lEta) = k2Data2%Beta(iOffE+iEta)
      EInv(lEta) = k2Data2%ZInv(iOffE+iEta)
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
  xpre(ij+1:ij+lZeta) = rKA(1:lZeta)*rKC(iEta)*sqrt(One/(Zeta(1:lZeta)+Eta(iEta)))
  ij = ij+lZeta
end do
if (ldot) then
  jPAO = 0
  do iPAO=1,mPAO
    PAO(jPAO+1:jPAO+lZeta*lEta) = xpre(1:lZeta*lEta)*PAO(jPAO+1:jPAO+lZeta*lEta)
    jPAO = jPAO+lZeta*lEta
  end do
end if
#ifdef _DEBUGPRINT_
call RecPrt(' PAO',' ',PAO,lZeta*lEta,mPAO)
#endif

return

end subroutine Screen_mck
