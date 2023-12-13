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
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Screen_g(iOffZ,iOffE,PAO,Scrtch,mPAO,nZeta,nEta,mZeta,mEta,lZeta,lEta, &
                    k2Data1,k2Data2, &
                    Zeta,ZInv,P,xA,xB, &
                    Eta,EInv,Q,xG,xD, &
                    iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2,CutGrd,l2DI, &
                    PreScr,nScrtch,IsChi,ChiI2)
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
!             April '92 modified for gradient estimate                 *
!***********************************************************************

use k2_structure, only: k2_type
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iOffZ, iOffE, mPAO, nZeta, nEta, mZeta, mEta, iphX1, iphY1, iphZ1, iphX2, iphY2, iphZ2, nScrtch, &
                                 IsChi
real(kind=wp), intent(inout) :: PAO(mZeta*mEta*mPAO)
real(kind=wp), intent(out) :: Scrtch(nScrtch), Zeta(nZeta), ZInv(nZeta), P(nZeta,3), xA(nZeta), xB(nZeta), Eta(nEta), EInv(nEta), &
                              Q(nEta,3), xG(nEta), xD(nEta)
integer(kind=iwp), intent(out) :: lZeta, lEta
type(k2_type), intent(in) :: k2Data1, k2Data2
real(kind=wp), intent(in) :: CutGrd, ChiI2
logical(kind=iwp), intent(in) :: l2DI, PreScr
integer(kind=iwp) :: i, iab, iabcd, icd, iDMin, iEP, iEta, ij, iMin, iOff, ip, ip1, ip2, iPAO, ipE, ipFac, ipOAP, ipP, ipPAO, &
                     ipZ, iZE, iZeta, jPAO, jPZ, l1, l2, nab, ncd
real(kind=wp) :: alpha, beta, Cut2, eMin, Et, Px, Py, Pz, qEta, Qx, Qy, Qz, qZeta, rEta, rKAB, rKCD, rqEta, rqZeta, rZeta, temp, &
                 vMax, zMin, Zt
logical(kind=iwp) :: ZPreScr, EPreScr
real(kind=wp), pointer :: ab(:,:), abG(:,:), cd(:,:), cdG(:,:)
real(kind=wp), external :: DNrm2_

nab = size(k2Data1%abG,1)/nZeta
ncd = size(k2Data2%abG,1)/nEta

ab(1:nZeta,1:nab) => k2Data1%abG(1:nZeta*nab,1)
abG(1:nZeta,1:nab) => k2Data1%abG(1:nZeta*nab,2)
cd(1:nEta,1:ncd) => k2Data2%abG(1:nEta*ncd,1)
cdG(1:nEta,1:ncd) => k2Data2%abG(1:nEta*ncd,2)

#ifdef _DEBUGPRINT_
if (l2DI) then
  call RecPrt(' ab   ',' ',ab,nZeta,nab)
  call RecPrt(' cd   ',' ',cd,nEta,ncd)
  call RecPrt(' abg  ',' ',abg,nZeta,nab)
  call RecPrt(' cdg  ',' ',cdg,nEta,ncd)
end if
call RecPrt('2nd order density matrix',' ',PAO,mZeta*mEta,mPAO)
#endif
if (PreScr .and. (.not. l2DI)) then
  write(u6,*) ' Screen: .not.l2DI no activated  prescr=',prescr,'  l2di=',l2di
  call Abend()
end if

Cut2 = CutGrd
lZeta = 0
lEta = 0
ip = 1

! Compute the prefactor. There are
! four sections of code here which will take care of the
! different cases of triangularization.

ipFac = ip
ip = ip+mZeta*mEta
ij = ipFac-1
do iEta=1,mEta
  Et = k2Data2%Zeta(iOffE+iEta)
  rKCD = k2Data2%Kappa(iOffE+iEta)
  do iZeta=1,mZeta
    Zt = k2Data1%Zeta(iOffZ+iZeta)
    rKAB = k2Data1%Kappa(iOffZ+iZeta)
    ij = ij+1
    if (IsChi == 1) then
      Scrtch(ij) = rKAB*rKCD*sqrt(One/(Zt+Et+(Zt*Et*ChiI2)*real(IsChi,kind=wp)))
    else
      Scrtch(ij) = rKAB*rKCD*sqrt(One/(Zt+Et))
    end if
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt(' Collected prefactors',' ',Scrtch(ipFac),mZeta,mEta)
#endif

! Modify the 2nd order density matrix with the prefactor.

ipPAO = ip
ip = ip+mZeta*mEta*mPAO
ipOAP = ip
ip = ip+mZeta*mEta*mPAO
jPAO = 1
do iPAO=1,mPAO
  do iZE=0,mZeta*mEta-1
    PAO(jPAO+iZE) = Scrtch(ipFac+iZE)*PAO(jPAO+iZE)
  end do
  jPAO = jPAO+mZeta*mEta
end do
#ifdef _DEBUGPRINT_
call RecPrt(' Modified 2nd order density matrix',' ',PAO,mZeta*mEta,mPAO)
#endif

! Scan the modified 2nd order density matrix for the
! largest absolute value.

vMax = DNrm2_(mZeta*mEta*mPAO,PAO,1)
! Skip prescreening if too low.
if (PreScr .and. (abs(vMax) < 0.5e-4_wp*CutGrd)) then
  lZeta = 0
  lEta = 0
  return
end if

! Assemble gradient estimate
! Multiply modified density matrix with the gradient estimate
! Compress all indices except zeta in ipZ
! Compress all indices except eta in ipE

if (PreScr) then

  ipZ = ip
  ip = ip+mZeta
  ipE = ip
  ip = ip+mEta
  do i=1,mZeta+mEta
    Scrtch(ipZ+i-1) = Zero
  end do
  do icd=1,ncd
    do iEta=1,mEta
      alpha = cd(iOffE+iEta,icd)
      beta = cdg(iOffE+iEta,icd)
      iZE = (iEta-1)*mZeta
      do iab=1,nab
        iabcd = (icd-1)*nab+iab
        iOff = (iabcd-1)*mZeta*mEta+iZE
        do iZeta=1,mZeta
          temp = (alpha*abg(iOffZ+iZeta,iab)+beta*ab(iOffZ+iZeta,iab))*PAO(iZeta+iOff)
          Scrtch(ipZ+iZeta-1) = Scrtch(ipZ+iZeta-1)+abs(temp)
          Scrtch(ipE+iEta-1) = Scrtch(ipE+iEta-1)+abs(temp)
        end do
      end do
    end do
  end do
else
  ipZ = 1
  ipE = 1
end if

nullify(ab,abG,cd,cdG)

if (ip-1 > nScrtch) then
  write(u6,*) ' Screen: ip-1 > nScrtch'
  write(u6,*) 'ip-1=',ip-1
  write(u6,*) 'nScrtch=',nScrtch
  call Abend()
end if

rEta = real(nEta*mPAO,kind=wp)
qEta = real(mEta*mPAO,kind=wp)
rqEta = rEta/qEta
ZPreScr = PreScr
if (PreScr) then
  iMin = iDMin(mZeta,Scrtch(ipZ),1)
  zMin = Scrtch(ipZ-1+iMin)
  if (zMin >= Cut2/rqEta) then
    ZPreScr = .false.
  end if
# ifdef _DEBUGPRINT_
  call RecPrt(' Screening array(Eta)',' ',Scrtch(ipZ),mZeta,1)
# endif
end if

rZeta = real(nZeta*mPAO,kind=wp)
qZeta = real(mZeta*mPAO,kind=wp)
rqZeta = rZeta/qZeta
EPreScr = PreScr
if (PreScr) then
  iMin = iDMin(mEta,Scrtch(ipE),1)
  eMin = abs(Scrtch(ipE-1+iMin))
  if (eMin >= Cut2/rqZeta) then
    EPreScr = .false.
  end if
# ifdef _DEBUGPRINT_
  call RecPrt(' Screening array(Zeta)',' ',Scrtch(ipE),mEta,1)
# endif
end if

! Prescreen Zeta

Px = One
Py = One
Pz = One
if (iphX1 /= 1) Px = -Px
if (iphY1 /= 1) Py = -Py
if (iphZ1 /= 1) Pz = -Pz
lZeta = 0
if (ZPreScr) then
  do iZeta=1,mZeta
    if (Scrtch(ipZ+iZeta-1) >= Cut2/rqEta) then
      lZeta = lZeta+1
      Zeta(lZeta) = k2Data1%Zeta(iOffZ+iZeta)
      P(lZeta,1) = k2Data1%PCoor(iOffZ+iZeta,1)*Px
      P(lZeta,2) = k2Data1%PCoor(iOffZ+iZeta,2)*Py
      P(lZeta,3) = k2Data1%PCoor(iOffZ+iZeta,3)*Pz
      xA(lZeta) = k2Data1%Alpha(iOffZ+iZeta)
      xB(lZeta) = k2Data1%Beta(iOffZ+iZeta)
      ZInv(lZeta) = k2Data1%ZInv(iOffZ+iZeta)
      ip2 = ipOAP+mEta*mPAO*(lZeta-1)-1
      do iEP=1,mEta*mPAO
        Scrtch(ip2+iEP) = PAO((iEP-1)*mZeta+iZeta)
      end do
    end if
  end do
else
  do iZeta=1,mZeta
    lZeta = lZeta+1
    Zeta(lZeta) = k2Data1%Zeta(iOffZ+iZeta)
    P(lZeta,1) = k2Data1%PCoor(iOffZ+iZeta,1)*Px
    P(lZeta,2) = k2Data1%PCoor(iOffZ+iZeta,2)*Py
    P(lZeta,3) = k2Data1%PCoor(iOffZ+iZeta,3)*Pz
    xA(lZeta) = k2Data1%Alpha(iOffZ+iZeta)
    xB(lZeta) = k2Data1%Beta(iOffZ+iZeta)
    ZInv(lZeta) = k2Data1%ZInv(iOffZ+iZeta)
  end do
  if (EPreScr) call DGeTMO(PAO,mZeta,mZeta,mEta*mPAO,Scrtch(ipOAP),mEta*mPAO)
end if
if (lZeta == 0) return

! Prescreen Eta

Qx = One
Qy = One
Qz = One
if (iphX2 /= 1) Qx = -Qx
if (iphY2 /= 1) Qy = -Qy
if (iphZ2 /= 1) Qz = -Qz
lEta = 0
if (EPreScr) then
  do iEta=1,mEta
    if (Scrtch(ipE+iEta-1) >= Cut2/rqZeta) then
      lEta = lEta+1
      Eta(lEta) = k2Data2%Zeta(iOffE+iEta)
      Q(lEta,1) = k2Data2%PCoor(iOffE+iEta,1)*Qx
      Q(lEta,2) = k2Data2%PCoor(iOffE+iEta,2)*Qy
      Q(lEta,3) = k2Data2%PCoor(iOffE+iEta,3)*Qz
      xG(lEta) = k2Data2%Alpha(iOffE+iEta)
      xD(lEta) = k2Data2%Beta(iOffE+iEta)
      EInv(lEta) = k2Data2%ZInv(iOffE+iEta)
      ip1 = ipOAP+iEta-1
      ip2 = ipPAO+(lEta-1)*mPAO*lZeta-1
      do jPZ=1,mPAO*lZeta
        Scrtch(ip2+jPZ) = Scrtch((jPZ-1)*mEta+ip1)
      end do
    end if
  end do
  ipP = ipPAO
  l1 = mPAO
  l2 = lZeta*lEta
else
  do iEta=1,mEta
    lEta = lEta+1
    Eta(lEta) = k2Data2%Zeta(iOffE+iEta)
    Q(lEta,1) = k2Data2%PCoor(iOffE+iEta,1)*Qx
    Q(lEta,2) = k2Data2%PCoor(iOffE+iEta,2)*Qy
    Q(lEta,3) = k2Data2%PCoor(iOffE+iEta,3)*Qz
    xG(lEta) = k2Data2%Alpha(iOffE+iEta)
    xD(lEta) = k2Data2%Beta(iOffE+iEta)
    EInv(lEta) = k2Data2%ZInv(iOffE+iEta)
  end do
  ipP = ipOAP
  l1 = lEta*mPAO
  l2 = lZeta
end if
if (lEta == 0) return

! Pick up the screened two-particle density.
! .not. PreScr : density is already in PAO

! Transpose mPAO,zeta,eta to zeta,eta,mPAO   or
!           eta,mPAO,zeta to zeta,eta,mPAO

if (PreScr) then
  if (ZPreScr .or. EPreScr) then
    call DGeTMO(Scrtch(ipP),l1,l1,l2,PAO,l2)
  end if
end if

#ifdef _DEBUGPRINT_
call RecPrt(' PAO',' ',PAO,lZeta*lEta,mPAO)
#endif

return

end subroutine Screen_g
