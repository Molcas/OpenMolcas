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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine NrmClc(Vec,lth,SubNam,MatNam)
!***********************************************************************
!                                                                      *
!     purpose: compute and print out norms for debuging purposes       *
!                                                                      *
!     input:                                                           *
!       Vec     : vector whose norm is going to be computed (lth)      *
!       SumNam  : name of a subroutine it is called from               *
!       MatNam  : name of a matrix                                     *
!       DeBug   : T - print out norm                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: lth
real(kind=wp), intent(in) :: Vec(lth)
character(len=*), intent(in) :: SubNam, MatNam
integer(kind=iwp) :: i
real(kind=wp) :: Q, R, S

R = sum(Vec(:)**2)
Q = sum(Vec(:))
S = Zero
do i=1,lth
  S = S+Vec(i)*real(i,kind=wp)
end do
write(u6,'(5A,3ES24.16,I8)') ' Norm of ',MatNam,' in ',SubNam,' = ',R,Q,S,lth

return

end subroutine NrmClc
