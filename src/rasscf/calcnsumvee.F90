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

use rasscf_global, only: lRoots
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: RotMat(lroots,lroots), DDG(lRoots,lRoots,lRoots,lRoots)
real(kind=wp) :: CalcNSumVee
real(kind=wp), allocatable :: Vee(:)

call mma_allocate(Vee,lRoots)
call CalcVee(Vee,RotMat,DDg)
CalcNSumVee = sum(Vee(:))
call mma_deallocate(Vee)

end function CalcNSumVee
