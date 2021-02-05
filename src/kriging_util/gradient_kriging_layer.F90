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

subroutine Gradient_Kriging_Layer(qInt,Grad,nInter)

implicit none
#include "stdalloc.fh"
integer nInter
real*8 qInt(nInter), Grad(nInter)
real*8, allocatable :: qInt_s(:), Grad_s(:)

call mma_allocate(qInt_s,nInter,Label='qInt_s')
call mma_allocate(Grad_s,nInter,Label='Grad_s')

call Trans_K(qInt,qInt_s,nInter,1)
call Gradient_Kriging(qInt_s,Grad_s,nInter)
call BackTrans_K(Grad_s,Grad,nInter,1)

call mma_deallocate(Grad_s)
call mma_deallocate(qInt_s)

end subroutine Gradient_Kriging_Layer
