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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

subroutine Energy_Kriging_Layer(qInt,Energy,nInter)

implicit none
#include "stdalloc.fh"
integer nInter
real*8 qInt(nInter), Energy
real*8, allocatable :: qInt_s(:)

call mma_allocate(qInt_s,nInter,Label='qInt_s')

call Trans_K(qInt,qInt_s,nInter,1)
call Energy_Kriging(qInt_s,Energy,nInter)

call mma_deallocate(qInt_s)

end subroutine Energy_Kriging_Layer
