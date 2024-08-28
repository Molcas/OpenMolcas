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

use Basis_Info, only: nCnttp, DBSC
use Center_Info, only: DC
use Symmetry_Info, only: nIrrep
use Constants, only: Zero
use Disp, only: IndDsp
use Definitions, only: wp

implicit none
#include "Molcas.fh"
#include "SysDef.fh"
integer nGrad
real*8 CGrad(3,MxAtom)
real*8 GradIn(nGrad)
character(len=LENIN5) CNames(MxAtom)
integer iCen
real*8, parameter :: tol = 1.0e-8_wp
logical, external :: TF
integer mdc, iIrrep, iCnttp, iCnt, iCo, kOp, nDisps, iCar, iComp
integer, external :: iPrmt, NrOpr
real*8 Xr

mdc = 0
iIrrep = 0

call dcopy_(3*MxAtom,[Zero],0,CGrad,1)
iCen = 0
!nCnttp_Valence = 0
!do iCnttp=1,nCnttp
!  if (dbsc(iCnttp)%Aux) Go To 999
!  nCnttp_Valence = nCnttp_Valence+1
!end do
!999 continue

do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      do iCo=0,nIrrep/dc(mdc)%nStab-1
        kop = dc(mdc)%iCoSet(iCo,0)
        nDispS = IndDsp(mdc,iIrrep)
        iCen = iCen+1
        do iCar=0,2
          iComp = 2**iCar
          if (TF(mdc,iIrrep,iComp)) then
            nDispS = nDispS+1
            XR = real(iPrmt(NrOpr(kop),icomp),kind=wp)
            CGrad(iCar+1,iCen) = XR*GradIn(nDispS)
          end if
        end do
        CNames(iCen) = dc(mdc)%LblCnt
      end do
    end do
  end if
end do

end subroutine TrGrd_Alaska
