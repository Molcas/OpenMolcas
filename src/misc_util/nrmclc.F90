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

implicit real*8(a-h,o-z)
real*8 Vec(lth)
character SubNam*(*), MatNam*(*)

R = DDot_(lth,Vec,1,Vec,1)
Q = DDot_(lth,[1.0d0],0,Vec,1)
S = 0.0d0
do i=1,lth
  S = S+Vec(i)*dble(i)
end do
write(6,'(5A,3E24.16,I8)') ' Norm of ',MatNam,' in ',SubNam,' = ',R,Q,S,lth

return

end subroutine NrmClc
