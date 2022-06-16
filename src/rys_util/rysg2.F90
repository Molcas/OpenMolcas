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
! Copyright (C) 1990,1991, Roland Lindh                                *
!               1990, IBM                                              *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Rysg2(iAnga,nRys,nT,Alpha,Beta,Gmma,Delta,Zeta,ZInv,nZeta,Eta,EInv,nEta,P,lP,Q,lQ,Coori,Coora,CoorAC,Array,nArray, &
                 Tvalue,ModU2,Cff2D,PAO,nPAO,Hess,nHess,IfGrd,IndGrd,IfHss,IndHss,nOp,iuvwx,IfG,mVec,Index_Out,lGrad,lHess,Tr)
!***********************************************************************
!                                                                      *
! Object: to compute the gradient of the two-electron integrals.       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified to 1st order derivatives October '91            *
!             Modified to 2nd order derviatives Mars '95 by            *
!             Anders Bernhardsson Theoretical Chemistry,               *
!             University of Lund                                       *
!***********************************************************************
!   @parameter iAnga     Angular momenta for each center
!   @parameter nRys      Order of Rys polynomia
!   @parameter nT        Number of alpha-beta-gamma-delta multiplies
!   @parameter Alpha     Exponents on 1st center
!   @parameter Beta      Exponents on 2nd center
!   @parameter Gmma      Exponents on 3rd center
!   @parameter Delta     Exponents on 4th center
!   @parameter Zeta      Alpha*Beta
!   @parameter Zeta      Zeta inverse
!   @parameter nZeta     Alpha Beta multiplies
!   @parameter Eta       Gmma*Delta
!   @parameter Eta       Eta inverse
!   @parameter nEta      Gmma Delta multiplies
!   @parameter P
!   @parameter lP        Length of P
!   @parameter Q
!   @parameter lQ        Length of Q
!   @parameter Coori     Coordinates of center just used to check AeqB CeqD etc.
!   @parameter Coora     Coordinates of center used in hor. recursion
!   @parameter CoorAC    Coordinates of center <max(la,lb),max(lc,ld)>
!   @parameter Array     Scratch and output area for 1st derivatives
!   @parameter nArray    Size of scratch area
!   @parameter PAO       Density
!   @parameter nPAO      Length of density
!   @parameter Hess      Output area for Hessian (added)
!   @parameter nHess     Size of Hessian
!   @parameter IfGrad    True for all 1st derivatives that are needed
!   @parameter IndGrad   Index in gradient on which integrals should be added
!   @parameter IfHss     True for all 2nd derivatives that are needed
!   @parameter IndHss    Index in Hess on which contracted integrals should be added
!   @parameter nOp       Operator number for the operator that generates center
!   @parameter iuvwx     Number of stabilizers
!   @parameter IfG       True for all centers on which derivatives should be calculated
!   @parameter Index_Out Index where first derivatives are stored (out)
!   @parameter lHess     True if 2nd derivatives should be calculated
!   @parameter lGrad     True if 1st derivatives should be calculated
!   @parameter Tr        True for all centers on which should be calculated via translation invariance

use vRys_RW, only: nMxRys
use Symmetry_Info, only: nIrrep, iOper
use Gateway_Info, only: ChiI2
use Gateway_global, only: IsChi, NoTab
use Index_Functions, only: nTri_Elem1
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), nRys, nT, nZeta, nEta, lP, lQ, nArray, nPAO, nHess, IndGrd(3,4,0:7), nOp(4), iuvwx(4)
real(kind=wp), intent(in) :: Alpha(nZeta), Beta(nZeta), Gmma(nEta), Delta(nEta), Zeta(nZeta), ZInv(nZeta), Eta(nEta), EInv(nEta), &
                             P(lP,3), Q(lQ,3), Coori(3,4), Coora(3,4), CoorAC(3,2), PAO(nT,nPAO)
real(kind=wp), intent(inout) :: Array(nArray), Hess(nHess)
external :: Tvalue, ModU2, Cff2D
logical(kind=iwp), intent(inout) :: IfGrd(3,4), IfHss(4,3,4,3), IfG(4), Tr(4)
integer(kind=iwp), intent(inout) :: IndHss(4,3,4,3,0:7)
integer(kind=iwp), intent(out) :: mVec, Index_Out(3,4)
logical(kind=iwp), intent(in) :: lGrad, lHess
integer(kind=iwp) :: i, iCent, iEta, Index1(3,4), Index2(3,4,4), Index3(3,3), Index4(2,6,3), iOff, ip, ip2D0, ip2D1, ip2D2, ipB00, &
                     ipB01, ipB10, ipDiv, ipEInv, ipEta, ipg2, ipP, ipPAQP, ipQ, ipQCPQ, ipScr, ipScr2, ipTmp, ipTv, ipU2, ipWgh, &
                     ipZeta, ipZInv, iStop, iZeta, jCar, JndGrd(3,4,0:7), kCent, la, lab, labMax, lb, lB00, lB01, lB10, lc, lCar, &
                     lcd, ld, lla, llb, llc, lld, lOp(4), MemFinal, n2D0, n2D1, n2D2, nabMax, ncdMax, ng(3), nh(3), nTR
logical(kind=iwp) :: KfGrd(3,4)
external :: Exp_1, Exp_2

KfGrd(:,:) = .false.
lOp(1) = iOper(nOp(1))
lOp(2) = iOper(nOp(2))
lOp(3) = iOper(nOp(3))
lOp(4) = iOper(nOp(4))
la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
lla = 0
llb = 0
llc = 0
lld = 0
if (any(IfGrd(:,1))) lla = 1
if (any(IfGrd(:,2))) llb = 1
if (any(IfGrd(:,3))) llc = 1
if (any(IfGrd(:,4))) lld = 1
do i=1,3
  if (IfHss(1,i,1,i)) lla = 2
  if (IfHss(2,i,2,i)) llb = 2
  if (IfHss(3,i,3,i)) llc = 2
  if (IfHss(4,i,4,i)) lld = 2
end do
lab = max(lla,llb)
lcd = max(llc,lld)
nabMax = la+lb+lab
ncdMax = lc+ld+lcd

! Allocate memory for the integral gradients.

ip = 1

MemFinal = 9*nT*nTri_Elem1(la)*nTri_Elem1(lb)*nTri_Elem1(lc)*nTri_Elem1(ld)
ip = ip+MemFinal

! Allocate memory for the 2D-integrals.

ip2D0 = ip
n2D0 = max((nabMax+1)*(ncdMax+1),(la+3)*(lb+3)*(ncdMax+1),(la+3)*(lb+3)*(lc+3)*(ld+3))
ip = ip+n2D0*3*nT*nRys

! Allocate memory for the 2nd order derivatives of the 2D-integrals.

ip2D1 = ip
n2D1 = max((nabMax+1)*(ncdMax+1),(la+3)*(lb+3)*(ncdMax+1),(la+1)*(lb+1)*(lc+1)*(ld+1)*3)
ip = ip+n2D1*3*nT*nRys

! Allocate memory for the coefficients in the recurrence relation
! of the 2D-integrals.

nTR = nT*nRys
ipPAQP = ip
ip = ip+nTR*3
ipQCPQ = ip
ip = ip+nTR*3
ipB10 = ip
lB10 = max(min(nabMax-1,1),0)
ip = ip+nTR*3*lB10
labMax = min(nabMax,ncdMax)
ipB00 = ip
lB00 = max(min(labMax,1),0)
ip = ip+nTR*3*lB00
ipB01 = ip
lB01 = max(min(ncdMax-1,1),0)
ip = ip+nTR*3*lB01
! Allocate memory for the roots.
ipU2 = ip
ip = ip+nT*nRys
! Allocate memory for Zeta, ZInv, Eta, EInv
ipZeta = ip
ip = ip+nT
ipEta = ip
ip = ip+nT
ipZInv = ip
ip = ip+nT
ipEInv = ip
ip = ip+nT
! Allocate memory for P and Q
ipP = ip
ip = ip+3*nT
ipQ = ip
ip = ip+3*nT
! Allocate memory for the inverse.
ipDiv = ip
ip = ip+nT
! Allocate memory for the arguments.
ipTv = ip
ip = ip+nT
if (ip-1 > nArray) then
  call WarningMessage(2,'Rysg2: ip-1 > nArray (pos. 1)')
  write(u6,*) 'ip,nArray=',ip,nArray
  call Abend()
end if

! Expand Zeta, ZInv, Eta, EInv, P, and Q

do iEta=1,nEta
  iOff = (iEta-1)*nZeta
  call dcopy_(nZeta,Zeta,1,Array(iOff+ipZeta),1)
  call dcopy_(nZeta,ZInv,1,Array(iOff+ipZInv),1)
  call dcopy_(nZeta,P(1,1),1,Array(iOff+ipP),1)
  iOff = iOff+nZeta*nEta
  call dcopy_(nZeta,P(1,2),1,Array(iOff+ipP),1)
  iOff = iOff+nZeta*nEta
  call dcopy_(nZeta,P(1,3),1,Array(iOff+ipP),1)
end do
do iZeta=1,nZeta
  iOff = iZeta-1
  call dcopy_(nEta,Eta,1,Array(iOff+ipEta),nZeta)
  call dcopy_(nEta,EInv,1,Array(iOff+ipEInv),nZeta)
  call dcopy_(nEta,Q(1,1),1,Array(iOff+ipQ),nZeta)
  iOff = iOff+nZeta*nEta
  call dcopy_(nEta,Q(1,2),1,Array(iOff+ipQ),nZeta)
  iOff = iOff+nZeta*nEta
  call dcopy_(nEta,Q(1,3),1,Array(iOff+ipQ),nZeta)
end do

! Compute the arguments for which we will compute the roots and the weights.

call Tvalue(Array(ipZeta),Array(ipEta),Array(ipP),Array(ipQ),nT,Array(ipTv),Array(ipDiv),IsChi,ChiI2)

! Compute roots and weights. Make sure that the weights ends up in
! the array where the z component of the 2D integrals will be.
! Call vRysRW if roots and weights are tabulated in various Taylor
! expansions. If not tabulated call RtsWgh.

! Pointer to z-component of 2D-integrals where the weights will be
! put directly. This corresponds to xyz2D(1,1,3,0,0).
ipWgh = ip2D0+2*nT*nRys
if ((nRys > nMxRys) .or. NoTab) then
  if (ip-1 > nArray) then
    call WarningMessage(2,'Rysg2: ip-1 > nArray (pos. 2)')
    write(u6,*) 'ip,nArray=',ip,nArray
    call Abend()
  end if

  call RtsWgh(Array(ipTv),nT,Array(ipU2),Array(ipWgh),nRys)
else
  if (ip-1 > nArray) then
    call WarningMessage(2,'Rysg2: ip-1 > nArray (pos. 3)')
    write(u6,*) 'ip,nArray=',ip,nArray
    call Abend()
  end if

  ! Make sure rys11/her11 is called  (la+1)
  call vRysRW(la+1,lb,lc,ld,Array(ipTv),Array(ipU2),Array(ipWgh),nT,nRys)
end if
! Drop ipTv
ip = ip-nT

! Modify the roots.

call ModU2(Array(ipU2),nT,nRys,Array(ipDiv))

! Drop ipDiv
ip = ip-nT

! Compute coefficients for the recurrence relations of the 2D-integrals

call Cff2D(max(nabMax-1,0),max(ncdMax-1,0),nRys,Array(ipZeta),Array(ipZInv),Array(ipEta),Array(ipEInv),nT,Coori,CoorAC,Array(ipP), &
           Array(ipQ),la+lab,lb,lc+lcd,ld,Array(ipU2),Array(ipPAQP),Array(ipQCPQ),Array(ipB10),Array(ipB00),labMax,Array(ipB01))
! Drop ipU2
ip = ip-nT*nRys
! Let go of Zeta, ZInv, Eta, and EInv
ip = ip-nT*4
! Let go of P and Q
ip = ip-6*nT

! -------------------------------
! Compute the intermediate 2D-integrals from the roots and weights.

call Rs2Dmm(Array(ip2D0),nT,nRys,nabMax,ncdMax,Array(ipPAQP),Array(ipQCPQ),Array(ipB10),Array(ipB00),Array(ipB01),la,lb,lc,ld, &
            IfHss,IfGrd)

! Drop ipB01
ip = ip-nTR*3*lB01
! Drop ipB00
ip = ip-nTR*3*lB00
! Drop ipB10
ip = ip-nTR*3*lB10
! Drop ipQCPQ
ip = ip-nTR*3
! Drop ipPAQP
ip = ip-nTR*3

! Apply the transfer equation to the intermediate 2D-integrals.

call HrrCtl_mck(Array(ip2D0),n2D0,Array(ip2D1),n2D1,la,lb,lc,ld,nabmax,ncdmax,nTR,Coora(1,1),Coora(1,2),Coora(1,3),Coora(1,4), &
                IfHss,IfGrd)

! Compute the gradients of the 2D-integrals. Copy some information
! which will be modified. This has to be done in order to facilitate
! partitioning.

ip2D2 = ip
n2D2 = (la+1)*(lb+1)*(lc+1)*(ld+1)*18
ip = ip+n2D2*nT*nRys
ipScr = ip
ip = ip+nT*nRys
ipTmp = ip
ip = ip+nT*nRys
ipScr2 = ip
ip = ip+nT*nRys
if (ip-1 > nArray) then
  call WarningMessage(2,'Rysg2: ip-1 > nArray (pos. 4)')
  write(u6,*) 'ip,nArray=',ip,nArray
  call Abend()
end if
JndGrd(:,:,:) = IndGrd(:,:,:)
do iCent=1,4
  do jCar=1,3
    do kCent=1,i
      if (iCent == kCent) then
        iStop = 3
      else
        iStop = jCar-1
      end if
      do lCar=1,istop
        if (jCar /= lCar) then
          if (IfHss(iCent,jCar,kCent,lCar)) then
            KfGrd(jCar,iCent) = .true.
            KfGrd(lCar,kCent) = .true.
          end if
        end if
      end do
    end do
  end do
end do
KfGrd(:,:) = KfGrd .or. IfGrd

call Rs2Dgh(Array(ip2D0),nT,nRys,la,lb,lc,ld,Array(ip2D1),Array(ip2D2),IfHss,IndHss,KfGrd,JndGrd,IfG,Coora,Alpha,Beta,Gmma,Delta, &
            nZeta,nEta,Array(ipScr),Array(ipScr2),Array(ipTmp),Index1,Index2,Index3,Index4,ng,nh,Exp_1,Exp_2,nZeta,nEta,nIrrep,Tr)
! Drop ipScr
ip = ip-nTR
! Drop ipScr2
ip = ip-nTR
! Drop ipTmp
ip = ip-nTR

! G2

ipg2 = ip
ip = ip+78
if (ip-1 > nArray) then
  call WarningMessage(2,'Rysg2: ip-1 > nArray (pos. 5)')
  write(u6,*) 'ip,nArray=',ip,nArray
  call Abend()
end if

if (lGrad) then
  KfGrd(:,:) = KfGrd .and. IfGrd
  call Assg1_mck(Array,nT,nRys,la,lb,lc,ld,Array(ip2D0),Array(ip2D1),KfGrd,Index1,mVec,Index_out)

end if

if (lHess) then
  call Assg2(Array(ipg2),nT,nRys,la,lb,lc,ld,Array(ip2D0),Array(ip2D1),Array(ip2D2),IfHss,Index3,Index4,ng,nh,PAO)

  ! Distribute the contributions to the molecular static hessian

  call Distg2(Array(ipg2),Hess,nHess,JndGrd,IfHss,IndHss,iuvwx,lOp,nop,Tr,IfG)
end if

! Drop ipg2
ip = ip-78
! Drop ip2D1
ip = ip-nTR*3*n2D1
! Drop ip2D0
ip = ip-nTR*3*n2D0
! Drop ip2D2
ip = ip-n2D2*nT*nRys

if (ip /= 1+MemFinal) then
  call WarningMessage(2,'Rysg2: ip /= 1+MemFinal (pos. 5)')
  write(u6,*) 'ip,MemFinal=',ip,MemFinal
  call Abend()
end if
IfGrd(:,:) = KfGrd

return

end subroutine Rysg2
