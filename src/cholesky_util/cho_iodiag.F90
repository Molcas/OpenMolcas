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

subroutine CHO_IODIAG(DIAG,IOPT)

implicit none
real*8 DIAG(*)
integer IOPT
character*7 FNAME
parameter(FNAME='CHODIAG')

call CHO_IODIAG_1(DIAG,IOPT,FNAME)

end subroutine CHO_IODIAG
