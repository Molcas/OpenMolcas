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
!***********************************************************************

subroutine Rysg1(iAnga,nRys,nT,Alpha,Beta,Gmma,Delta, &
                 Zeta,ZInv,nZeta,Eta,EInv,nEta, &
                 P,lP,Q,lQ,Coori,Coora,CoorAC, &
                 Array,nArray,Tvalue,ModU2,Cff2D,PAO,nPAO,Grad,nGrad,IfGrad,IndGrd,kOp,iuvwx)
!***********************************************************************
!                                                                      *
! Object: to compute the gradient of the two-electron integrals.       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified to 1st order derivatives October '91            *
!***********************************************************************

use vRys_RW, only: nMxRys
use Symmetry_Info, only: iOper
use Gateway_Info, only: ChiI2
use Gateway_global, only: IsChi, NoTab
use Breit, only: nOrdOp
use Definitions, only: wp, iwp
#if defined (_DEBUGPRINT_) || defined (_CHECK_)
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), nRys, nT, nZeta, nEta, lP, lQ, nArray, nPAO, nGrad, IndGrd(3,4), kOp(4), iuvwx(4)
real(kind=wp), intent(in) :: Alpha(nZeta), Beta(nZeta), Gmma(nEta), Delta(nEta), Zeta(nZeta), ZInv(nZeta), Eta(nEta), EInv(nEta), &
                             P(lP,3), Q(lQ,3), Coori(3,4), Coora(3,4), CoorAC(3,2), PAO(nT,nPAO)
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Array(nArray)
external :: Tvalue, ModU2, Cff2D
logical(kind=iwp), intent(in) :: IfGrad(3,4)
integer(kind=iwp) :: iEta, Indx(3,4), iOff, ip, ip2D0, ip2D1, ipB00, ipB01, ipB10, ipDiv, ipEInv, ipEta, ipP, ipPAQP, ipQ, ipQCPQ, &
                     ipScr, ipTmp, ipTv, ipU2, ipWgh, ipZeta, ipZInv, iZeta, JndGrd(3,4), la, lab, labMax, lb, lB00, lB01, lB10, &
                     lc, lcd, ld, lla, llb, llc, lld, lOp(4), mVec, n2D0, n2D1, nabMax, ncdMax, nTR
real(kind=wp) :: Temp(9)
logical(kind=iwp) :: JfGrad(3,4)
external :: Exp_1, Exp_2

lOp(1) = iOper(kOp(1))
lOp(2) = iOper(kOp(2))
lOp(3) = iOper(kOp(3))
lOp(4) = iOper(kOp(4))
#ifdef _DEBUGPRINT_
call RecPrt(' In Rysg1:P',' ',P(1:nZeta,:),nZeta,3)
call RecPrt(' In Rysg1:Q',' ',Q(1:nEta,:),nEta,3)
call RecPrt(' In Rysg1:Zeta',' ',Zeta,nZeta,1)
call RecPrt(' In Rysg1:ZInv',' ',ZInv,nZeta,1)
call RecPrt(' In Rysg1:Eta',' ',Eta,nEta,1)
call RecPrt(' In Rysg1:EInv',' ',EInv,nEta,1)
call RecPrt(' In Rysg1:Alpha',' ',Alpha,nZeta,1)
call RecPrt(' In Rysg1:Beta ',' ',Beta,nZeta,1)
call RecPrt(' In Rysg1:Gamma',' ',Gmma,nEta,1)
call RecPrt(' In Rysg1:Delta',' ',Delta,nEta,1)
call RecPrt(' In Rysg1:Coora',' ',Coora,3,4)
call RecPrt(' In Rysg1:Coori',' ',Coori,3,4)
call RecPrt(' In Rysg1:CoorAC',' ',CoorAC,3,2)
write(u6,*) ' In Rysg1: iAnga=',iAnga
#endif
la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
lla = 0
llb = 0
llc = 0
lld = 0
if (any(IfGrad(:,1))) lla = 1
if (any(IfGrad(:,2))) llb = 1
if (any(IfGrad(:,3))) llc = 1
if (any(IfGrad(:,4))) lld = 1
lab = max(lla,llb)
lcd = max(llc,lld)
nabMax = la+lb+lab
ncdMax = lc+ld+lcd

! Allocate memory for the integral gradients.

ip = 1
!ipAC = ip
!ip = ip+nT*nPAO*9

! Allocate memory for the 2D-integrals.

ip2D0 = ip
n2D0 = max((nabMax+1)*(ncdMax+1),(la+2)*(lb+2)*(ncdMax+1),(la+2)*(lb+2)*(lc+2)*(ld+2))
ip = ip+n2D0*3*nT*nRys

! Allocate memory for the 1st order derivatives of the 2D-integrals.

ip2D1 = ip
n2D1 = max((nabMax+1)*(ncdMax+1),(la+2)*(lb+2)*(ncdMax+1),(la+1)*(lb+1)*(lc+1)*(ld+1)*3)
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
!#define _CHECK_
#ifdef _CHECK_
if (ip-1 > nArray) then
  call WarningMessage(2,'Rysg1: ip-1 =/= nArray (pos.1)')
  write(u6,*) ' nArray=',nArray
  write(u6,*) ' ip-1  =',ip-1
  write(u6,*) ' nRys  =',nRys
  write(u6,*) ' nZeta =',nZeta
  write(u6,*) ' nEta  =',nEta
  call Abend()
end if
#endif

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

call Tvalue(Array(ipZeta),Array(ipEta),Array(ipP),Array(ipQ),nT,Array(ipTv),Array(ipDiv),IsChi,ChiI2,nOrdOp)

! Compute roots and weights. Make sure that the weights ends up in
! the array where the z component of the 2D integrals will be.
! Call vRysRW if roots and weights are tabulated in various Taylor
! expansions. If not tabulated call RtsWgh.

! Pointer to z-component of 2D-integrals where the weights will be
! put directly. This corresponds to xyz2D(1,1,3,0,0).
ipWgh = ip2D0+2*nT*nRys
if ((nRys > nMxRys) .or. NoTab) then
# ifdef _CHECK_
  if (ip-1 > nArray) then
    call WarningMessage(2,'Rysg1: ip-1 =/= nArray (pos.2)')
    write(u6,*) ' nArray=',nArray
    write(u6,*) ' ip-1  =',ip-1
    call Abend()
  end if
# endif
  call RtsWgh(Array(ipTv),nT,Array(ipU2),Array(ipWgh),nRys,nOrdOp)
else
# ifdef _CHECK_
  if (ip-1 > nArray) then
    call WarningMessage(2,'Rysg1: ip-1 =/= nArray (pos.3)')
    write(u6,*) ' nArray=',nArray
    write(u6,*) ' ip-1  =',ip-1
    call Abend()
  end if
# endif
  call vRysRW(la+1,lb,lc,ld,Array(ipTv),Array(ipU2),Array(ipWgh),nT,nRys,nOrdOp)
end if
! Drop ipTv
ip = ip-nT

! Modify the roots.

call ModU2(Array(ipU2),nT,nRys,Array(ipDiv))
! Drop ipDiv
ip = ip-nT

! Compute coefficients for the recurrence relations of the 2D-integrals

call Cff2D(max(nabMax-1,0),max(ncdMax-1,0),nRys,Array(ipZeta),Array(ipZInv),Array(ipEta),Array(ipEInv),nT,Coori,CoorAC,Array(ipP), &
           Array(ipQ),la+lab,lb,lc+lcd,ld,Array(ipU2),Array(ipPAQP),Array(ipQCPQ),Array(ipB10),Array(ipB00),labMax,Array(ipB01), &
           nOrdOp)
! Drop ipU2
ip = ip-nT*nRys
! Let go of Zeta, ZInv, Eta, and EInv
ip = ip-nT*4
! Let go of P and Q
ip = ip-6*nT

! Compute the intermediate 2D-integrals from the roots and weights.

call vRys2Dm(Array(ip2D0),nT,nRys,nabMax,ncdMax,Array(ipPAQP),Array(ipQCPQ),Array(ipB10),Array(ipB00),Array(ipB01),la,lb,lc,ld, &
             IfGrad)
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

call HrrCtl(Array(ip2D0),n2D0,Array(ip2D1),n2D1,la,lb,lc,ld,nabmax,ncdmax,nTR,Coora(1,1),Coora(1,2),Coora(1,3),Coora(1,4),IfGrad)

! Compute the gradients of the 2D-integrals. Copy some information
! which will be modified. This has to be done in order to facilitate
! partitioning.

ipScr = ip
ip = ip+nT*nRys
ipTmp = ip
ip = ip+nT
JndGrd(:,:) = IndGrd
JfGrad(:,:) = IfGrad
call Rys2Dg(Array(ip2D0),nT,nRys,la,lb,lc,ld,Array(ip2D1),JfGrad,JndGrd,Coora,Alpha,Beta,Gmma,Delta,nZeta,nEta,Array(ipScr), &
            Array(ipTmp),Indx,Exp_1,Exp_2,nZeta,nEta)
! Drop ipScr
ip = ip-nTR
! Drop ipTmp
ip = ip-nT

! Assemble the gradients of the ERI's

call Assg1(Temp,PAO,nT,nRys,la,lb,lc,ld,Array(ip2D0),Array(ip2D1),JfGrad,Indx,mVec)
! Drop ip2D1
ip = ip-nTR*3*n2D1
! Drop ip2D0
ip = ip-nTR*3*n2D0

! Distribute the contributions to the molecular gradient

call Distg1(Temp,Grad,nGrad,JfGrad,JndGrd,iuvwx,lOp)
! Drop ipAC
!ip = ip-nT*nPAO*9
#ifdef _CHECK_
if (ip /= 1) then
  call WarningMessage(2,'Rysg1: ip=/=1')
  call Abend()
end if
#endif

return

end subroutine Rysg1
