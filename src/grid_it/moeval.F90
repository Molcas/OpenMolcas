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
!               2000, Valera Veryazov                                  *
!               2014, Thomas Dresselhaus                               *
!***********************************************************************

subroutine MOEval(MOValue,nMOs,nCoor,CCoor,CMOs,nCMO,DoIt,nDrv,mAO)
!***********************************************************************
!      Author:Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN. November 1995                           *
!                                                                      *
!      Modified: Thomas Dresselhaus, March 2014                        *
!                Added ability to calculate 2nd derivative as well     *
!***********************************************************************

use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Phase_Info, only: iPhase
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
! nDrv: Between 0 and 2. The highest derivative to be calc.
! mAO:  Memory slots per point and basis functions. Should be >=1 for nDrv=0, >=4 for nDrv=1, >=10 for nDrv=2.
integer(kind=iwp), intent(in) :: nMOs, nCoor, nCMO, DoIt(nMOs), nDrv, mAO
real(kind=wp), intent(out) :: MOValue(mAO,nCoor,nMOs)
real(kind=wp), intent(in) :: CCoor(3,nCoor), CMOs(nCMO)
#include "WrkSpc.fh"
#include "print.fh"
integer(kind=iwp) :: iAng, iAO, iAOttp, iBas, iCmp, iCnt, iCnttp, iDrv, iG, ipAng, ipAOs, ipRadial, iPrim, iPrint, ipSOs, iptmp, &
                     ipx, ipxyz, ipy, ipz, iRout, iScrt1, iScrt2, iShll, iSkal, kSh, mdc, mRad, nAngular, nAO, nCnt, nDeg, nForm, &
                     nOp, nRadial, nScr1, nScr2, nSO, nTerm, nTest, ntmp, nxyz
real(kind=wp) :: A(3), px, py, pz, RA(3), Thr
integer(kind=iwp), external :: NrOpr

! Statement functions
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

iRout = 112
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(u6,*) ' In MOEval'
end if
MOValue(:,:,:) = Zero

! Loop over shells.

iSkal = 0
Thr = Zero

do iAng=S%iAngMx,0,-1

  if (S%MaxPrm(iAng) == 0) goto 100
  if (S%MaxBas(iAng) == 0) goto 100

  ! Scratch area for contraction step

  nScr1 = S%MaxPrm(iAng)*nElem(iAng)
  call GetMem('Scrtch','ALLO','REAL',iScrt1,nScr1)

  ! Scratch area for the transformation to spherical gaussians

  nScr2 = S%MaxPrm(iAng)*nElem(iAng)
  call GetMem('ScrSph','Allo','Real',iScrt2,nScr2)

  ! Loop over basis sets. Skip if basis set do not include
  ! angular momentum functions as specified above.

  iAOttp = 0
  mdc = 0
  do iCnttp=1,nCnttp

    nTest = dbsc(iCnttp)%nVal
    if (iAng+1 > nTest) Go To 101
    if (dbsc(iCnttp)%Aux) Go To 101
    if (dbsc(iCnttp)%Frag) Go To 101
    nCnt = dbsc(iCnttp)%nCntr
    iShll = dbsc(iCnttp)%iVal+iAng
    iPrim = Shells(iShll)%nExp
    if (iPrim == 0) Go To 101
    iBas = Shells(iShll)%nBasis
    if (iBas == 0) Go To 101
    if (Shells(iShll)%Prjct) then
      iCmp = 2*iAng+1
    else
      iCmp = nElem(iAng)
    end if

    call OrdExpD2C(iPrim,Shells(iShll)%Exp,iBas,Shells(iShll)%pCff)
    kSh = dbsc(iCnttp)%iVal+iAng

    ! Loop over unique centers of basis set "iCnttp"

    do iCnt=1,nCnt

      iAO = iAOttp+(iCnt-1)*dbsc(iCnttp)%lOffAO+Shells(kSh)%kOffAO
      A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

      !------- Allocate memory for SO and AO values

      mRad = nDrv+1

      nForm = 0
      do iDrv=0,nDrv
        nForm = nForm+nElem(iDrv)
      end do
      nTerm = 2**nDrv

      nAO = (iCmp*iBas*nCoor)*(mAO)
      nSO = nAO*nIrrep/dc(mdc+iCnt)%nStab
      nDeg = nIrrep/dc(mdc+iCnt)%nStab
      call GetMem('AOs','Allo','Real',ipAOs,nAO)
      call GetMem('SOs','Allo','Real',ipSOs,nSO)
      call dcopy_(nSO,[Zero],0,Work(ipSOs),1)
      nxyz = nCoor*3*(iAng+mRad)
      call GetMem('xyz','Allo','Real',ipxyz,nxyz)
      ntmp = nCoor
      call GetMem('tmp','Allo','Real',iptmp,ntmp)
      nRadial = iBas*nCoor*mRad
      call GetMem('Radial','Allo','Real',ipRadial,nRadial)

      nAngular = 5*nForm*nTerm
      call GetMem('Angular','Allo','Inte',ipAng,nAngular)

      !------- Loops over symmetry operations operating on the basis set center.

      do iG=0,nIrrep/dc(mdc+iCnt)%nStab-1
        iSkal = iSkal+1
        call OA(dc(mdc+iCnt)%iCoSet(iG,0),A,RA)
        ipx = iPhase(1,dc(mdc+iCnt)%iCoSet(iG,0))
        ipy = iPhase(2,dc(mdc+iCnt)%iCoSet(iG,0))
        ipz = iPhase(3,dc(mdc+iCnt)%iCoSet(iG,0))
        px = real(iPhase(1,dc(mdc+iCnt)%iCoSet(iG,0)),kind=wp)
        py = real(iPhase(2,dc(mdc+iCnt)%iCoSet(iG,0)),kind=wp)
        pz = real(iPhase(3,dc(mdc+iCnt)%iCoSet(iG,0)),kind=wp)
        nOp = NrOpr(dc(mdc+iCnt)%iCoSet(iG,0))

        !----- Evaluate AOs at RA

        call dcopy_(nAO,[Zero],0,Work(ipAOs),1)
        call AOEval(iAng,nCoor,CCoor,Work(ipxyz),RA,Shells(iShll)%Transf,RSph(ipSph(iAng)),nElem(iAng),iCmp,iWork(ipAng),nTerm, &
                    nForm,Thr,mRad,iPrim,iPrim,Shells(iShll)%Exp,Work(ipRadial),iBas,Shells(iShll)%pCff,Work(ipAOs),mAO,px,py,pz, &
                    ipx,ipy,ipz)

        !----- Distribute contributions to the SOs

        call SOAdpt(Work(ipAOs),mAO,nCoor,iBas,iCmp,nOp,Work(ipSOs),nDeg,iAO)

      end do ! iG

      !------- Distribute contributions to the MOs

      call SODist(Work(ipSOs),mAO,nCoor,iBas,iCmp,nDeg,MOValue,nMOs,iAO,CMOs,nCMO,DoIt)

      call GetMem('Radial','Free','Real',ipRadial,nRadial)
      call GetMem('Angular','Free','Inte',ipAng,nAngular)
      call GetMem('tmp','Free','Real',iptmp,ntmp)
      call GetMem('xyz','Free','Real',ipxyz,nxyz)
      call GetMem('AOs','Free','Real',ipAOs,nAO)
      call GetMem('SOs','Free','Real',ipSOs,nSO)

    end do ! iCnt
101 continue
    mdc = mdc+dbsc(iCnttp)%nCntr
    iAOttp = iAOttp+dbsc(iCnttp)%lOffAO*dbsc(iCnttp)%nCntr
  end do
  call GetMem('ScrSph','Free','Real',iScrt2,nScr2)
  call GetMem('Scrtch','Free','Real',iScrt1,nScr1)
100 continue
end do

return

end subroutine MOEval
