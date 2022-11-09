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
! Copyright (C) 1987, Jeppe Olsen                                      *
!               1989, Markus P. Fuelscher                              *
!***********************************************************************

subroutine IWRTMA(IMAT,NROW,NCOL,MAXROW,MAXCOL)
! AUTHOR: J. OLSEN, UNIV. OF LUND, SWEDEN, APRIL 1987
! MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
!                M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
!
! WRITE A MATRIX OF INTEGER VALUES

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MAXROW, MAXCOL, IMAT(MAXROW,MAXCOL), NROW, NCOL
integer(kind=iwp) :: I

do I=1,NROW
  write(u6,1110) IMAT(I,1:NCOL)
end do

return

1110 format(/,1X,8I10,/,(1X,8I10))

end subroutine IWRTMA
