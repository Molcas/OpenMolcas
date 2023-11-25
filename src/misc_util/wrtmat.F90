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

subroutine WRTMAT(A,NROW,NCOL,NMROW,NMCOL)
! AUTHOR: J. OLSEN, UNIV. OF LUND, SWEDEN, APRIL 1987
! MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
!                M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
!
! PRINT REAL VALUED MATRIX

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NROW, NCOL, NMROW, NMCOL
real(kind=wp), intent(in) :: A(NMROW,NMCOL)
integer(kind=iwp) :: I

do I=1,NROW
  write(u6,1010) I,A(I,1:NCOL)
end do

return

1010 format('0',I3,2X,4(ES15.8),/,(6X,4(ES15.8)))

end subroutine WRTMAT
