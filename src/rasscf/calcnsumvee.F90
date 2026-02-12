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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

function CalcNSumVee(RotMat,DDg)

use stdalloc, only: mma_allocate, mma_deallocate
use rasscf_global, only: lRoots

implicit none
real*8, dimension(lRoots,lRoots,lRoots,lRoots) :: DDG
real*8, dimension(lroots,lroots) :: RotMat
real*8, dimension(:), allocatable :: Vee
real*8 CalcNSumVee
integer IState
#include "warnings.h"

call mma_allocate(Vee,lRoots)
CalcNSumVee = 0.0d0
call CalcVee(Vee,RotMat,DDg)
do IState=1,lRoots
  CalcNSumVee = CalcNSumVee+Vee(IState)
end do
call mma_deallocate(Vee)

end function CalcNSumVee
