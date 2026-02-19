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

subroutine CMSMatRot(Mat,A,I,J,N)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: I, J, N
real(kind=wp) :: Mat(N,N), A
real(kind=wp), allocatable :: TM(:,:)

call mma_allocate(TM,2,N,Label='TM')
TM(1,:) = Mat(I,:)
TM(2,:) = Mat(J,:)
Mat(J,:) = cos(A)*TM(2,:)+sin(A)*TM(1,:)
Mat(I,:) = -sin(A)*TM(2,:)+cos(A)*TM(1,:)
call mma_deallocate(TM)

end subroutine CMSMatRot
