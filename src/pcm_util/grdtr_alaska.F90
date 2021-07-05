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
!               2002, Roland Lindh                                     *
!***********************************************************************

subroutine GrdTr_Alaska(GradIn,MxAto,GradOut,nGrad)
!*******************************************************************
!                                                                  *
!  The inverse of                                                  *
!  Transforms a symmetry adapted gradient to unsymmetric  form     *
!                                                                  *
!   Written by Anders Bernhardsson                                 *
!   960427                                                         *
!   Modified by Roland Lindh                                       *
!   020115                                                         *
!                                                                  *
!*******************************************************************

use Basis_Info
use Center_Info
use Symmetry_Info, only: nIrrep
implicit real*8(a-h,o-z)
parameter(tol=1d-8)
#include "Molcas.fh"
#include "disp.fh"
#include "real.fh"
#include "WrkSpc.fh"
real*8 GradIn(3,MxAto), GradOut(nGrad)
logical, external :: TF

iIrrep = 0

mdc = 0
iCen = 1
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux) return
  do iCnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    nDispS = IndDsp(mdc,iIrrep)

    do iCar=0,2
      iComp = 2**iCar
      if (TF(mdc,iIrrep,iComp)) then
        nDispS = nDispS+1
        GradOut(nDispS) = GradIn(iCar+1,iCen)
      end if
    end do
    iCen = iCen+nIrrep/dc(mdc)%nStab
  end do
end do

return

end subroutine GrdTr_Alaska
