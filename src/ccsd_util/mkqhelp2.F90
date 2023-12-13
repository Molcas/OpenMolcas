!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine mkqhelp2(vector,dimv,length,factor)
! this routine does vector = vector*factor
! vector - multiplied vector (I/O)
! dimv   - dimension of vector
! length - length of vector to be multiplied
! factor - scaling factor
!
! $N.B. this routine should be substituted by mv0s3v

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimv, length
real(kind=wp), intent(inout) :: vector(dimv)
real(kind=wp), intent(in) :: factor

if (length > 0) vector(1:length) = vector(1:length)*factor

return

end subroutine mkqhelp2
