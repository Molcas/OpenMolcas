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
! Copyright (C) 1995,2001, Roland Lindh                                *
!***********************************************************************

subroutine PCMgrd1(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,final,nZeta,la,lb,A,RB,nRys,Array,nArr,Ccoor,nOrdOp,Grad,nGrad, &
                   IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,iStabM,nStabM)
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, May '95                                 *
!                                                                      *
!             Modified to PCM gradients September 2001, Lund, by       *
!             R. Lindh.                                                *
!***********************************************************************

use PCM_arrays, only: PCMTess
use Center_Info

implicit real*8(A-H,O-Z)
external TNAI1, Fake, XCff2D
#include "real.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "rctfld.fh"
integer IndGrd(3,2), kOp(2), lOper(nComp), iStabM(0:nStabM-1), iDCRT(0:7)
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6), Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta), rKappa(nZeta), &
       P(nZeta,3), A(3), RB(3), CCoor(3,nComp), Array(nZeta*nArr), Grad(nGrad), DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2)
logical IfGrad(3,2)
! Local arrys

real*8 C(3), TC(3), Coori(3,4), CoorAC(3,2)
logical NoLoop, JfGrad(3,4)
integer iAnga(4), iStb(0:7), JndGrd(3,4), lOp(4), iuvwx(4)
character ChOper(0:7)*3
data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

iRout = 151
iPrint = nPrint(iRout)

! Modify the density matrix with the prefactor

nDAO = nElem(la)*nElem(lb)
do iDAO=1,nDAO
  do iZeta=1,nZeta
    Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
    DAO(iZeta,iDAO) = Fact*DAO(iZeta,iDAO)
  end do
end do
if (iPrint >= 99) call RecPrt('DAO',' ',DAO,nZeta,nDAO)

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
ipDAO = nip
nip = nip+nAlpha*nBeta*nElem(la)*nElem(lb)*nElem(nOrdOp)
if (nip-1 > nZeta*nArr) then
  write(6,*) 'nip-1 > nZeta*nArr'
  call ErrTra()
  call AbEnd()
end if
nArray = nZeta*nArr-nip+1

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = nOrdOp
iAnga(4) = 0
call dcopy_(3,A,1,Coori(1,1),1)
call dcopy_(3,RB,1,Coori(1,2),1)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
lOp(1) = kOp(1)
lOp(2) = kOp(2)

ipAOff = ipA
do iBeta=1,nBeta
  call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
  ipAOff = ipAOff+nAlpha
end do

ipBOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
  ipBOff = ipBOff+1
end do

! Loop over the tiles

! pcm_solvent remove the loop
!do iTs=1,nTs
do iTs=1,1
! pcm_solvent end
  ! pcm_solvent put "charge" to 1
  !Q = PCM_SQ(1,iTs)+PCM_SQ(2,iTs)
  Q = One
  ! pcm_solvent end
  NoLoop = Q == Zero
  if (NoLoop) Go To 111
  ! Pick up the tile coordinates
  C(1:3) = PCMTess(1:3,iTs)

  if (iPrint >= 99) call RecPrt('C',' ',C,1,3)

  ! Generate stabilizer of C

  nStb = 1
  iStb(0) = 0

  ! Find the DCR for M and S

  call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
  Fact = -dble(nStabM)/dble(LmbdT)

  if (iPrint >= 99) then
    write(6,*) ' Q=',Q
    write(6,*) ' Fact=',Fact
    call RecPrt('DAO*Fact*Q',' ',Array(ipDAO),nZeta*nDAO,nElem(nOrdOp))
    write(6,*) ' m      =',nStabM
    write(6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
    write(6,*) ' s      =',nStb
    write(6,'(9A)') '(S)=',(ChOper(iStb(ii)),ii=0,nStb-1)
    write(6,*) ' LambdaT=',LmbdT
    write(6,*) ' t      =',nDCRT
    write(6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
  end if
  iuvwx(3) = nStb
  iuvwx(4) = nStb
  call ICopy(6,IndGrd,1,JndGrd,1)
  do i=1,3
    do j=1,2
      JfGrad(i,j) = IfGrad(i,j)
    end do
  end do

  ! No derivatives with respect to the third or fourth center.
  ! The positions of the points in the external field are frozen.

  call ICopy(3,[0],0,JndGrd(1,3),1)
  JfGrad(1,3) = .false.
  JfGrad(2,3) = .false.
  JfGrad(3,3) = .false.
  call ICopy(3,[0],0,JndGrd(1,4),1)
  JfGrad(1,4) = .false.
  JfGrad(2,4) = .false.
  JfGrad(3,4) = .false.
  mGrad = 0
  do iCar=1,3
    do i=1,2
      if (JfGrad(iCar,i)) mGrad = mGrad+1
    end do
  end do
  if (iPrint >= 99) write(6,*) ' mGrad=',mGrad
  if (mGrad == 0) Go To 111

  do lDCRT=0,nDCRT-1
    lOp(3) = NrOpr(iDCRT(lDCRT))
    lOp(4) = lOp(3)
    call OA(iDCRT(lDCRT),C,TC)
    call dcopy_(3,TC,1,CoorAC(1,2),1)
    call dcopy_(3,TC,1,Coori(1,3),1)
    call dcopy_(3,TC,1,Coori(1,4),1)

    call DYaX(nZeta*nDAO,Fact*Q,DAO,1,Array(ipDAO),1)

    ! Compute integrals with the Rys quadrature.

    nT = nZeta
    nDiff = 1
    mRys = (la+lb+2+nDiff+nOrdOp)/2
    call Rysg1(iAnga,mRys,nT,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,Coori,Coori,CoorAC, &
               Array(nip),nArray,TNAI1,Fake,XCff2D,Array(ipDAO),nDAO*nElem(nOrdOp),Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)

    !call RecPrt(' In PCMgrd:Grad',' ',Grad,nGrad,1)
  end do  ! End loop over DCRs

111 continue
end do     ! End loop over centers in the external field

!call GetMem(' Exit PCMgrd','LIST','REAL',iDum,iDum)
return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(final)
  call Unused_integer(nRys)
  call Unused_real_array(Ccoor)
  call Unused_integer_array(lOper)
end if

end subroutine PCMgrd1
