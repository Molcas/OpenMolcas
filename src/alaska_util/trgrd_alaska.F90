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
! Copyright (C) 1996, Anders Bernhardsson                              *
!               2006, Giovanni Ghigo                                   *
!***********************************************************************

subroutine TrGrd_Alaska(CGrad,CNames,GradIn,nGrad,iCen)
!***********************************************************************
!                                                                      *
!      Transforms a symmetry adapted gradient to unsymmetric form.     *
!                                                                      *
!       Written by Anders Bernhardsson                                 *
!       Adapted by Giovanni Ghigo                                      *
!       University of Torino, July 2006                                *
!***********************************************************************

use Basis_Info, only: DBSC, nCnttp
use Center_Info, only: DC
use Symmetry_Info, only: nIrrep
use Disp, only: IndDsp
use Molcas, only: LenIn, MxAtom
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(out) :: CGrad(3,MxAtom)
character(len=LenIn+5), intent(out) :: CNames(MxAtom)
real(kind=wp), intent(in) :: GradIn(nGrad)
integer(kind=iwp), intent(out) :: iCen
integer(kind=iwp) :: iCar, iCnt, iCnttp, iCo, iComp, kOp, mdc, nDisps
integer(kind=iwp), parameter :: iIrrep = 0
real(kind=wp), parameter :: tol = 1.0e-8_wp
integer(kind=iwp), external :: iPrmt, NrOpr
logical(kind=iwp), external :: TF

CGrad(:,:) = Zero
iCen = 0
!nCnttp_Valence = 0
!do iCnttp=1,nCnttp
!  if (dbsc(iCnttp)%Aux) exit
!  nCnttp_Valence = nCnttp_Valence+1
!end do

mdc = 0
do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      do iCo=0,nIrrep/dc(mdc)%nStab-1
        kOp = dc(mdc)%iCoSet(iCo,0)
        nDispS = IndDsp(mdc,iIrrep)
        iCen = iCen+1
        do iCar=0,2
          iComp = 2**iCar
          if (TF(mdc,iIrrep,iComp)) then
            nDispS = nDispS+1
            CGrad(iCar+1,iCen) = real(iPrmt(NrOpr(kOp),iComp),kind=wp)*GradIn(nDispS)
          end if
        end do
        CNames(iCen) = dc(mdc)%LblCnt
      end do
    end do
  end if
end do

end subroutine TrGrd_Alaska
