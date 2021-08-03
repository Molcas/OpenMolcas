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
!               2003, Michael Diedenhofen                              *
!***********************************************************************

subroutine COSGrd( &
#                 define _CALLING_
#                 include "grd_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of electronic COSMO cont. *
!         integrals.                                                   *
!                                                                      *
!             M. Diedenhofen Nov. 2003                                 *
!             changes pcmgrd routines which do not take into account   *
!             the contribution of a non fixed grid                     *
!     Author of orig routine: Roland Lindh,                            *
!             Dept. of Theoretical Chemistry, University               *
!             of Lund, Sweden, May '95                                 *
!                                                                      *
!             Modified to PCM gradients September 2001, Lund, by       *
!             R. Lindh.                                                *
!***********************************************************************

use PCM_arrays, only: PCM_SQ, PCMTess
use Center_Info, only: dc
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
integer(kind=iwp) :: i, iAlpha, iAnga(4), iBeta, iCar, iDAO, iDCRT(0:7), ii, iIrrep, ipA, ipAOff, ipB, ipBOff, ipDAO, iPrint, &
                     iRout, iStb(0:7), iTs, iuvwx(4), iZeta, j, JndGrd(3,4), kat, lDCRT, LmbdT, lOp(4), mGrad, mRys, nArray, nDAO, &
                     nDCRT, nDisp, nip, nRys, nStb, nT
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), Fact, Q, TC(3)
logical(kind=iwp) :: JfGrad(3,4)
character(len=3), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: NrOpr
external :: Cff2D, Fake, TNAI1
#include "Molcas.fh"
#include "print.fh"
#include "disp.fh"
#include "rctfld.fh"
! Statement function for Cartesian index
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

#include "macros.fh"
unused_var(rFinal)
unused_var(Ccoor)
unused_var(lOper)

iRout = 151
iPrint = nPrint(iRout)

nRys = nHer

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
  call WarningMessage(2,'Error in COSGrd')
  write(u6,*) 'nip-1 > nZeta*nArr'
  call Abend()
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

do iTs=1,nTs
  Q = PCM_SQ(1,iTs)
  if (Q == Zero) cycle
  ! Pick up the tile coordinates
  C(1) = PCMTess(1,iTs)
  C(2) = PCMTess(2,iTs)
  C(3) = PCMTess(3,iTs)
  kat = nint(PCMTess(4,iTs))
  if (iPrint >= 99) call RecPrt('C',' ',C,3,1)

  ! Generate stabilizor of C

  nStb = 1
  iStb(0) = 0

  ! Find the DCR for M and S

  call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
  Fact = -real(nStabM,kind=wp)/real(LmbdT,kind=wp)

  if (iPrint >= 99) then
    write(u6,*) ' Q=',Q
    write(u6,*) ' Fact=',Fact
    call RecPrt('DAO*Fact*Q',' ',Array(ipDAO),nZeta*nDAO,nElem(nOrdOp))
    write(u6,*) ' m      =',nStabM
    write(u6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
    write(u6,*) ' s      =',nStb
    write(u6,'(9A)') '(S)=',(ChOper(iStb(ii)),ii=0,nStb-1)
    write(u6,*) ' LambdaT=',LmbdT
    write(u6,*) ' t      =',nDCRT
    write(u6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
  end if
  iuvwx(3) = nStb
  iuvwx(4) = nStb

  call ICopy(6,IndGrd,1,JndGrd,1)
  do i=1,3
    do j=1,2
      JfGrad(i,j) = IfGrad(i,j)
    end do
  end do
  ! skip contributions if the segment belongs to the atom
  if (mdc == kat) then
    ! skip 1 center
    JfGrad(1,1) = .false.
    JfGrad(2,1) = .false.
    JfGrad(3,1) = .false.
  end if
  if (ndc == kat) then
    ! skip 2 center
    JfGrad(1,2) = .false.
    JfGrad(2,2) = .false.
    JfGrad(3,2) = .false.
  end if

  ! up third center

  iIrrep = 0
  nDisp = IndDsp(kat,iIrrep)
  do iCar=0,2
    nDisp = nDisp+1
    JndGrd(iCar+1,1) = abs(JndGrd(iCar+1,1))
    JndGrd(iCar+1,2) = abs(JndGrd(iCar+1,2))
    JndGrd(iCar+1,3) = -nDisp
    JfGrad(iCar+1,3) = .false.
  end do

  call ICopy(3,[0],0,JndGrd(1,4),1)
  JfGrad(1,4) = .false.
  JfGrad(2,4) = .false.
  JfGrad(3,4) = .false.

  mGrad = 0
  do iCar=1,3
    do i=1,3
      if (JfGrad(iCar,i)) mGrad = mGrad+1
    end do
  end do
  if (iPrint >= 99) write(u6,*) ' mGrad=',mGrad
  if (mGrad == 0) cycle

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
    mRys = nRys

    call Rysg1(iAnga,mRys,nT,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,Coori,Coori,CoorAC, &
               Array(nip),nArray,TNAI1,Fake,Cff2D,Array(ipDAO),nDAO*nElem(nOrdOp),Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)

    !call RecPrt(' In COSgrd:Grad',' ',Grad,nGrad,1)
  end do  ! End loop over DCRs

end do    ! End loop over centers in the external field

return

end subroutine COSGrd
