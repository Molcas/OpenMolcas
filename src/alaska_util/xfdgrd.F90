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
! Copyright (C) 1995, Roland Lindh                                     *
!***********************************************************************

subroutine XFdGrd( &
#                 define _CALLING_
#                 include "grd_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, May 1995                                *
!***********************************************************************

use external_centers, only: iXPolType, nOrd_XF, nXF, XF
use Center_Info, only: dc
use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
integer(kind=iwp) :: i, iAlpha, iAnga(4), iBeta, iCar, iChxyz, iDAO, iDCRT(0:7), iDum, iFd, ii, iOrdOp, ipA, ipAOff, ipB, ipBOff, &
                     ipDAO, iPrint, iRout, iStb(0:7), iuvwx(4), iZeta, j, jCoSet(8,8), JndGrd(3,4), jpDAO, lDCRT, LmbdT, lOp(4), &
                     mGrad, mRys, nArray, nDAO, nDCRT, nDiff, nip, nStb, nT
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), Fact, TC(3), TZFd(3), ZFd(3)
logical(kind=iwp) :: JfGrad(3,4), NoLoop
character(len=3), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: iChAtm, NrOpr
external :: Fake, TNAI1, XCff2D
#include "print.fh"

#include "macros.fh"
unused_var(rFinal)
unused_var(nHer)
unused_var(Ccoor(1))
unused_var(nOrdOp)
unused_var(nComp)

iRout = 151
iPrint = nPrint(iRout)

! Modify the density matrix with the prefactor

nDAO = nTri_Elem1(la)*nTri_Elem1(lb)
do iDAO=1,nDAO
  do iZeta=1,nZeta
    Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
    DAO(iZeta,iDAO) = Fact*DAO(iZeta,iDAO)
  end do
end do
if (iPrint >= 99) call RecPrt('DAO',' ',DAO,nZeta,nDAO)

! Loop over charges and dipole moments in the external field

if ((nOrd_XF > 1) .or. (iXPolType > 0)) then
  call WarningMessage(2,'Error in xfdgrd')
  write(u6,*) 'Sorry, gradients are not implemented for'
  write(u6,*) 'higher XF than dipoles or for polarisabilities'
  call Quit_OnUserError()
end if

do iOrdOp=0,nOrd_XF

  nip = 1
  ipA = nip
  nip = nip+nAlpha*nBeta
  ipB = nip
  nip = nip+nAlpha*nBeta
  ipDAO = nip
  nip = nip+nAlpha*nBeta*nTri_Elem1(la)*nTri_Elem1(lb)*nTri_Elem1(iOrdOp)
  if (nip-1 > nZeta*nArr) then
    call WarningMessage(2,'Error in xfdgrd')
    write(u6,*) 'nip-1 > nZeta*nArr'
    call Abend()
  end if
  nArray = nZeta*nArr-nip+1

  iAnga(1) = la
  iAnga(2) = lb
  iAnga(3) = iOrdOp
  iAnga(4) = 0
  Coori(:,1) = A(:)
  Coori(:,2) = RB(:)

  ! Find center to accumulate angular momentum on. (HRR)

  if (la >= lb) then
    CoorAC(:,1) = A(:)
  else
    CoorAC(:,1) = RB(:)
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

  ! Loop over centers of the external field.

  iDum = 0
  do iFd=1,nXF
    if (iOrdOp == 0) then
      ZFd(1) = XF(4,iFd)
      NoLoop = ZFd(1) == Zero
    else
      ZFd(1:3) = XF(5:7,iFd)
      NoLoop = (ZFd(1) == Zero) .and. (ZFd(2) == Zero) .and. (ZFd(3) == Zero)
    end if
    if (NoLoop) cycle
    ! Pick up the center coordinates
    C(1:3) = XF(1:3,iFd)

    if (iPrint >= 99) call RecPrt('C',' ',C,1,3)

    ! Generate stabilizer of C

    iChxyz = iChAtm(C)
    call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

    ! Find the DCR for M and S

    call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
    Fact = -real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    if (iPrint >= 99) then
      write(u6,*) ' ZFd=',(ZFd(i),i=1,nTri_Elem1(iOrdOp))
      write(u6,*) ' Fact=',Fact
      call RecPrt('DAO*Fact*ZFd()',' ',Array(ipDAO),nZeta*nDAO,nTri_Elem1(iOrdOp))
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
    JndGrd(:,1:2) = IndGrd(:,:)
    do i=1,3
      do j=1,2
        JfGrad(i,j) = IfGrad(i,j)
      end do
    end do

    ! No derivatives with respect to the third or fourth center.
    ! The positions of the points in the external field are frozen.

    JndGrd(:,3) = 0
    JfGrad(1,3) = .false.
    JfGrad(2,3) = .false.
    JfGrad(3,3) = .false.
    JndGrd(:,4) = 0
    JfGrad(1,4) = .false.
    JfGrad(2,4) = .false.
    JfGrad(3,4) = .false.
    mGrad = 0
    do iCar=1,3
      do i=1,2
        if (JfGrad(iCar,i)) mGrad = mGrad+1
      end do
    end do
    if (iPrint >= 99) write(u6,*) ' mGrad=',mGrad
    if (mGrad == 0) cycle

    do lDCRT=0,nDCRT-1
      lOp(3) = NrOpr(iDCRT(lDCRT))
      lOp(4) = lOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      CoorAC(:,2) = TC(:)
      Coori(:,3) = TC(:)
      Coori(:,4) = TC(:)

      if (iOrdOp == 0) then
        Array(ipDAO:ipDAO+nZeta*nDAO-1) = Fact*ZFd(1)*pack(DAO,.true.)
      else
        call OA(iDCRT(lDCRT),ZFd,TZFd)
        jpDAO = ipDAO
        Array(jpDAO:jpDAO+nZeta*nDAO-1) = Fact*TZFd(1)*pack(DAO,.true.)
        jpDAO = jpDAO+nZeta*nDAO
        Array(jpDAO:jpDAO+nZeta*nDAO-1) = Fact*TZFd(2)*pack(DAO,.true.)
        jpDAO = jpDAO+nZeta*nDAO
        Array(jpDAO:jpDAO+nZeta*nDAO-1) = Fact*TZFd(3)*pack(DAO,.true.)
      end if

      ! Compute integrals with the Rys quadrature.

      nT = nZeta
      nDiff = 1
      mRys = (la+lb+2+nDiff+iOrdOp)/2
      call Rysg1(iAnga,mRys,nT,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,Coori,Coori,CoorAC, &
                 Array(nip),nArray,TNAI1,Fake,XCff2D,Array(ipDAO),nDAO*nTri_Elem1(iOrdOp),Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)

      !call RecPrt(' In XFdGrd:Grad',' ',Grad,nGrad,1)
    end do ! End loop over DCRs

  end do   ! End loop over centers in the external field

end do     ! End loop over charges and dipole moments

return

end subroutine XFdGrd
