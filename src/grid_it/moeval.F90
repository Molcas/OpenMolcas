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
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
! nDrv: Between 0 and 2. The highest derivative to be calc.
! mAO:  Memory slots per point and basis functions. Should be >=1 for nDrv=0, >=4 for nDrv=1, >=10 for nDrv=2.
integer(kind=iwp), intent(in) :: nMOs, nCoor, nCMO, DoIt(nMOs), nDrv, mAO
real(kind=wp), intent(out) :: MOValue(mAO,nCoor,nMOs)
real(kind=wp), intent(in) :: CCoor(3,nCoor), CMOs(nCMO)
#include "print.fh"
integer(kind=iwp) :: iAng, iAO, iAOttp, iBas, iCmp, iCnt, iCnttp, iDrv, iG, iPrim, iPrint, ipx, ipy, ipz, iRout, iShll, iSkal, &
                     kSh, mdc, mRad, nAngular, nAO, nCnt, nDeg, nElem, nForm, nOp, nRadial, nSO, nTerm, nTest, nxyz
real(kind=wp) :: A(3), px, py, pz, RA(3), Thr
integer(kind=iwp), allocatable :: Ang(:)
real(kind=wp), allocatable :: AOs(:), Radial(:), SOs(:), xyz(:)
integer(kind=iwp), external :: NrOpr

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

  if (S%MaxPrm(iAng) == 0) cycle
  if (S%MaxBas(iAng) == 0) cycle

  ! Loop over basis sets. Skip if basis set do not include
  ! angular momentum functions as specified above.

  iAOttp = 0
  mdc = 0
  do iCnttp=1,nCnttp

    if (iCnttp > 1) then
      mdc = mdc+dbsc(iCnttp-1)%nCntr
      iAOttp = iAOttp+dbsc(iCnttp-1)%lOffAO*dbsc(iCnttp-1)%nCntr
    end if

    nTest = dbsc(iCnttp)%nVal
    if (iAng+1 > nTest) cycle
    if (dbsc(iCnttp)%Aux) cycle
    if (dbsc(iCnttp)%Frag) cycle
    nCnt = dbsc(iCnttp)%nCntr
    iShll = dbsc(iCnttp)%iVal+iAng
    iPrim = Shells(iShll)%nExp
    if (iPrim == 0) cycle
    iBas = Shells(iShll)%nBasis
    if (iBas == 0) cycle
    if (Shells(iShll)%Prjct) then
      iCmp = 2*iAng+1
    else
      iCmp = (iAng+1)*(iAng+2)/2
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
        nForm = nForm+(iDrv+1)*(iDrv+2)/2
      end do
      nTerm = 2**nDrv

      nAO = (iCmp*iBas*nCoor)*(mAO)
      nSO = nAO*nIrrep/dc(mdc+iCnt)%nStab
      nDeg = nIrrep/dc(mdc+iCnt)%nStab
      call mma_allocate(AOs,nAO,label='AOs')
      call mma_allocate(SOs,nSO,label='SOs')
      SOs(:) = Zero
      nxyz = nCoor*3*(iAng+mRad)
      call mma_allocate(xyz,nxyz,label='xyz')
      nRadial = iBas*nCoor*mRad
      call mma_allocate(Radial,nRadial,label='Radial')

      nAngular = 5*nForm*nTerm
      call mma_allocate(Ang,nAngular,label='Angular')

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

        AOs(:) = Zero
        nElem = (iAng+1)*(iAng+2)/2
        call AOEval(iAng,nCoor,CCoor,xyz,RA,Shells(iShll)%Transf,RSph(ipSph(iAng)),nElem,iCmp,Ang,nTerm,nForm,Thr,mRad,iPrim, &
                    iPrim,Shells(iShll)%Exp,Radial,iBas,Shells(iShll)%pCff,AOs,mAO,px,py,pz,ipx,ipy,ipz)

        !----- Distribute contributions to the SOs

        call SOAdpt(AOs,mAO,nCoor,iBas,iCmp,nOp,SOs,nDeg,iAO)

      end do ! iG

      !------- Distribute contributions to the MOs

      call SODist(SOs,mAO,nCoor,iBas,iCmp,nDeg,MOValue,nMOs,iAO,CMOs,nCMO,DoIt)

      call mma_deallocate(AOs)
      call mma_deallocate(SOs)
      call mma_deallocate(xyz)
      call mma_deallocate(Radial)
      call mma_deallocate(Ang)

    end do ! iCnt
  end do
end do

return

end subroutine MOEval
