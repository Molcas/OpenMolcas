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

implicit none
integer I, J, N
real*8 A
real*8, dimension(N,N) :: Mat, TM
integer K

do K=1,N
  TM(I,K) = Mat(I,K)
  TM(J,K) = Mat(J,K)
end do
do K=1,N
  Mat(J,K) = cos(A)*TM(J,K)+sin(A)*TM(I,K)
  Mat(I,K) = -sin(A)*TM(J,K)+cos(A)*TM(I,K)
end do

end subroutine CMSMatRot
