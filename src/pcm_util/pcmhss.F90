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
! Copyright (C) 1995,2001,2008, Roland Lindh                           *
!***********************************************************************

subroutine PCMHss( &
#                 define _CALLING_
#                 include "hss_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, May 1995                                *
!                                                                      *
!             Modified to PCM gradients September 2001, Lund, by       *
!             R. Lindh.                                                *
!             Modified to PCM Hessian February 2008, Lund by           *
!             R. Lindh.                                                *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use PCM_arrays, only: PCM_SQ, PCMTess
use Center_Info, only: dc
use rctfld_module, only: nTS
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "hss_interface.fh"
integer(kind=iwp) :: iAlpha, iAnga(4), iAtom, iBeta, iCar, iDAO, iDCRT(0:7), iIrrep, idx(3,4), ipA, ipAOff, ipB, ipBOff, ipDAO, &
                     iPrint, iRout, iStb(0:7), iTs, iuvwx(4), iZeta, jAtom, jCar, JndGrd(0:2,0:3,0:7), &
                     JndHss(0:3,0:2,0:3,0:2,0:7), lDCRT, LmbdT, mOp(4), mRys, nArray, nDAO, nDCRT, nDiff, nFinal, nip, nla, nlb, &
                     nOOp, nStb
real(kind=wp) :: Coori(3,4), CoorAC(3,2), C(3), EInv, Eta, Fact, TC(3), q_i
logical(kind=iwp) :: IfG(0:3), JfGrd(0:2,0:3), JfHss(0:3,0:2,0:3,0:2), NoLoop, Tr(0:3)
integer(kind=iwp), external :: NrOpr
external :: Fake, TNAI1, XCff2D
#include "print.fh"

#include "macros.fh"
unused_var(rFinal)
unused_var(nHer)
unused_var(Ccoor)
unused_var(lOper)

!                                                                      *
!***********************************************************************
!                                                                      *
! We will have five terms here!
!
! 1) Sum(i) q_i V_ie^xy
! 2) Sum(ij) V_in^y Q_ij V_je^x
! 3) Sum(ij) V_ie^y Q_ij V_je^x
! Maurizio to add comments for the last two terms!

iRout = 151
iPrint = nPrint(iRout)

nla = (la+1)*(la+2)/2
nlb = (lb+1)*(lb+2)/2
nOOp = (nOrdOp+1)*(nOrdOp+2)/2

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
ipDAO = nip
nip = nip+nAlpha*nBeta*nla*nlb*nOOp
if (nip-1 > nZeta*nArr) then
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
if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
mOp(1) = nOp(1)
mOp(2) = nOp(2)

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

! Modify the density matrix with the prefactor

nDAO = nla*nlb
do iDAO=1,nDAO
  do iZeta=1,nZeta
    Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
    DAO(iZeta,iDAO) = Fact*DAO(iZeta,iDAO)
  end do
end do
if (iPrint >= 99) call RecPrt('DAO',' ',DAO,nZeta,nDAO)

! Generate stabilizer of C, i.e. a center of a tile.

nStb = 1
iStb(0) = 0

! Loop over the tiles

do iTs=1,nTs
  q_i = PCM_SQ(1,iTs)+PCM_SQ(2,iTs)
  NoLoop = q_i == Zero
  if (NoLoop) cycle
  ! Pick up the tile coordinates
  C(1:3) = PCMTess(1:3,iTs)

  if (iPrint >= 99) call RecPrt('C',' ',C,1,3)
  call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
  Fact = -q_i*real(nStabM,kind=wp)/real(LmbdT,kind=wp)

  Array(ipDAO:ipDAO+nZeta*nDAO-1) = Fact*pack(DAO,.true.)

  iuvwx(3) = nStb
  iuvwx(4) = nStb

  do lDCRT=0,nDCRT-1
    mOp(3) = NrOpr(iDCRT(lDCRT))
    mOp(4) = mOp(3)
    call OA(iDCRT(lDCRT),C,TC)
    call dcopy_(3,TC,1,CoorAC(1,2),1)
    call dcopy_(3,TC,1,Coori(1,3),1)
    call dcopy_(3,TC,1,Coori(1,4),1)

    ! Initialize JfGrd, JndGrd, JfHss, and JndHss.

    JfGrd(:,:) = .false.
    JndGrd(:,:,0:nSym-1) = 0
    JfHss(:,:,:,:) = .false.
    JndHss(:,:,:,:,0:nSym-1) = 0

    ! Overwrite with information in IfGrd, IndGrd, IfHss,
    ! and IndHss. This sets up the info for the first two
    ! centers! Make sure that no translational invariance
    ! is used.

    do iAtom=0,1
      do iCar=0,2
        JfGrd(iCar,iAtom) = Ifgrd(iCar,iAtom)
        do iIrrep=0,nSym-1
          JndGrd(iCar,iAtom,iIrrep) = abs(IndGrd(iCar,iAtom,iIrrep))
        end do
        do jAtom=0,1
          do jCar=0,2
            JfHss(iAtom,iCar,jAtom,jCar) = IfHss(iAtom,iCar,jAtom,jCar)
            do iIrrep=0,nSym-1
              JndHss(iAtom,iCar,jAtom,jCar,iIrrep) = abs(IndHss(iAtom,iCar,jAtom,jCar,iIrrep))
            end do
          end do
        end do
      end do
    end do

    ! Derivatives with respect to the operator is computed via
    ! the translational invariance.
    ! Note: We want no such thing!
    !
    ! The third center is calculated by translational invariance.
    ! This requires the 2nd derivatives on the other centers.
    ! Note: We want no such thing!

    Tr(:) = .false.
    IfG(:) = [.true.,.true.,.false.,.false.]
    JfGrd(:,:) = .false.

    ! Compute integrals with the Rys quadrature.

    nDiff = 2
    mRys = (la+lb+2+nDiff+nOrdOp)/2
    Eta = One
    EInv = One
    nFinal = 0
    call Rysg2(iAnga,mRys,nZeta,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,Coori,Coori,CoorAC, &
               Array(nip),nArray,TNAI1,Fake,XCff2D,Array(ipDAO),nDAO*nOOp,Hess,nHess,JfGrd,JndGrd,JfHss,JndHss,mOp,iuvwx,IfG, &
               nFinal,idx,.false.,.true.,Tr)

    !call RecPrt(' In PCMHss:Hess',' ',Hess,nHess,1)
  end do  ! End loop over DCRs

end do    ! End loop over centers in the external field

return

end subroutine PCMHss
