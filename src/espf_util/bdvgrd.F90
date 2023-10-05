!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine BdVGrd( &
#                 define _CALLING_
#                 include "grd_interface.fh"
                 )

use Index_Functions, only: nTri_Elem1
use Center_Info, only: dc
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
integer(kind=iwp) :: i, iAlpha, iAnga(4), iBeta, iCar, iChxyz, iDAO, iDCRT(0:7), iDum, ii, iStb(0:7), iOrdOp, ipA, ipAOff, ipB, &
                     ipBOff, ipDAO, iPnt, IPotFl, iPrint, iuvwx(4), iZeta, jCoSet(8,8), JndGrd(3,4), jpDAO, lDCRT, LmbdT, lOp(4), &
                     mGrad, mRys, nArray, nDAO, nDCRT, nDiff, nGrdPt, nip, nStb, nT
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), Fact, TC(3), TZFd(3), ZFd(3), ZFdx, ZFdy, ZFdz
logical(kind=iwp) :: ESPFexist, JfGrad(3,4), NoLoop
character(len=180) :: Key
character(len=3), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: iChAtm, IsFreeUnit, NrOpr
character(len=180), external :: Get_Ln
external :: Fake, TNAI1, XCff2D
#include "macros.fh"
unused_var(rFinal(1,1,1,1,1))
unused_var(nHer)
unused_var(nOrdOp)
unused_var(nComp)

iPrint = 5

! Modify the density matrix with the prefactor

nDAO = nTri_Elem1(la)*nTri_Elem1(lb)
do iDAO=1,nDAO
  do iZeta=1,nZeta
    Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
    DAO(iZeta,iDAO) = Fact*DAO(iZeta,iDAO)
  end do
end do
if (iPrint >= 99) then
  call RecPrt('DAO',' ',DAO,nZeta,nDAO)
  call RecPrt('Alpha',' ',Alpha,1,nAlpha)
  call RecPrt('Beta ',' ',Beta,1,nBeta)
  call RecPrt('Zeta ',' ',Zeta,1,nZeta)
end if

! Loop over charges and dipole moments in the external field

iOrdOp = 0

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
ipDAO = nip
nip = nip+nAlpha*nBeta*nTri_Elem1(la)*nTri_Elem1(lb)*nTri_Elem1(iOrdOp)
if (nip-1 > nZeta*nArr) then
  write(u6,*) 'nip-1 > nZeta*nArr'
  call Abend()
end if
nArray = nZeta*nArr-nip+1

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = iOrdOp
iAnga(4) = 0
Coori(:,1) = A
Coori(:,2) = RB

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  CoorAC(:,1) = A
else
  CoorAC(:,1) = RB
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

! Loop over centers of the grid
! But how to retrieve the grid ???
! I just read it in the ESPF file !
call F_Inquire('ESPF.DATA',ESPFExist)
if (.not. ESPFExist) then
  write(u6,*) ' Error ! No ESPF.DATA file found'
  call Quit_OnUserError()
end if
IPotFl = IsFreeUnit(1)
call Molcas_Open(IPotFl,'ESPF.DATA')
do
  Key = Get_Ln(IPotFl)
  if (Key(1:10) == 'GRID      ') then
    call Get_I1(2,nGrdPt)
    exit
  end if
end do
close(IPotFl)

iDum = 0
do iPnt=1,nGrdPt
  ZFd(1) = CCoor((iPnt-1)*4+4)
  NoLoop = ZFd(1) == Zero
  if (NoLoop) cycle
  ! Pick up the center coordinates
  C(1) = CCoor((iPnt-1)*4+1)
  C(2) = CCoor((iPnt-1)*4+2)
  C(3) = CCoor((iPnt-1)*4+3)

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
  JndGrd(:,1:2) = IndGrd
  JfGrad(:,1:2) = IfGrad

  ! No derivatives with respect to the third or fourth center.
  ! The positions of the points in the external field are frozen.

  JndGrd(:,3:4) = 0
  JfGrad(:,3:4) = .false.
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
    CoorAC(:,2) = TC
    Coori(:,3) = TC
    Coori(:,4) = TC

    if (iOrdOp == 0) then
      call DYaX(nZeta*nDAO,Fact*ZFd(1),DAO,1,Array(ipDAO),1)
    else
      call OA(iDCRT(lDCRT),ZFd,TZFd)
      jpDAO = ipDAO
      ZFdx = TZFd(1)
      call DYaX(nZeta*nDAO,Fact*ZFdx,DAO,1,Array(jpDAO),1)
      jpDAO = jpDAO+nZeta*nDAO
      ZFdy = TZFd(2)
      call DYaX(nZeta*nDAO,Fact*ZFdy,DAO,1,Array(jpDAO),1)
      jpDAO = jpDAO+nZeta*nDAO
      ZFdz = TZFd(3)
      call DYaX(nZeta*nDAO,Fact*ZFdz,DAO,1,Array(jpDAO),1)
    end if

    ! Compute integrals with the Rys quadrature.

    nT = nZeta
    nDiff = 1
    mRys = (la+lb+2+nDiff+iOrdOp)/2
    call Rysg1(iAnga,mRys,nT,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,Coori,Coori,CoorAC, &
               Array(nip),nArray,TNAI1,Fake,XCff2D,Array(ipDAO),nDAO*nTri_Elem1(iOrdOp),Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)

    !call RecPrt(' In BdVGrd:Grad',' ',Grad,nGrad,1)
  end do  ! End loop over DCRs

end do  ! End loop over centers in the grid

return

end subroutine BdVGrd
